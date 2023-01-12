#!/usr/bin/env nextflow

///////////////////////////////
//// grab from samples.tsv ////
///////////////////////////////

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

fh = new File(params.inputFile)
def csv_content = fh.getText('utf-8')

def data_iterator = parseCsv(csv_content, separator: '\t', readFirstLine: false)

def fqArray = []
def f5Array = []
def rebasecallArray = []
def sampleArrayPairs = []
i = 0
//Sample-Codename FqPath  f5Path  ComparisonPair  ComparisonOrder
for (line in data_iterator) {
    fqArray[i] = [line[0],line[1]]
    f5Array[i] = [line[0],line[2]]
    rebasecallArray[i] = [line[3],line[0],line[1],line[2]]
    sampleArrayPairs[i] = [line[0],line[3],line[4]]
    i += 1;
}

fastqChannel = Channel.from(fqArray).unique()
signalChannel = Channel.from(f5Array).unique()
rebasecallChannel = Channel.from(rebasecallArray).unique()

workflow {
    // if rebase on same machine
    // reBasecall( rebasecallChannel )
    // reBasecall.out equiv to fastqChannel
    catFq( fastqChannel )
    miniMapping( params.genome, catFq.out )
    sam2bam( miniMapping.out )
    nanopolishEventalign( signalChannel.join(catFq.out).join(sam2bam.out) )
    xporeDataprep( nanopolishEventalign.out )
    ////// CHANNEL CHICANERY //////
    // basedata: ID, PAIR, AB
    basedata = Channel.from(sampleArrayPairs)
    // sampleAB.A: SAMPLE, NUMBER, A
    // sampleAB.B: SAMPLE, NUMBER, B 
    Channel
	.from(sampleArrayPairs)
	.branch {
		A: it[2] =~ 'A'
		B: it[2] =~ 'B'
	}
	.set { sampleAB }
    // SAMPLE, NUMBER, AB, PATH
    As = sampleAB.A.combine(xporeDataprep.out[0], by: 0)
    Bs = sampleAB.B.combine(xporeDataprep.out[0], by: 0)
    // pairChannel: NUMBER, SAMPLE, AB, PATH, SAMPLE, AB, PATH
    pairChannel = As.combine(Bs, by: 1)
    //////////////////////////////////
    makeYml( pairChannel )
    xporeDiffmod( makeYml.out )
    xporePostprocessing( xporeDiffmod.out[0] )
    // now for my train data prep //
    filterXpore( xporePostprocessing.out )
    // CK, XPOREPATH, AK, A, DPPATHA, BK, B, DPPATHB
    trainChan = filterXpore.out[0].combine(pairChannel, by: 0)
    xporeReadnames( trainChan )
    trainChan2 = xporeReadnames.out.combine(rebasecallChannel, by: 0)
    xporeTraindata( trainChan2 )
}

process reBasecall {
    input:
    tuple val(comparisonpair), val(runname), path(sample_fq), path(sample_f5)

    output:
    tuple val(runname), path(sample_fq)

    script:
    """
    if guppy_basecaller --version | grep -q "6.4"; then
        mkdir $sample_fq
        guppy_basecaller --flowcell FLO-PRO002 --kit SQK-RNA002 -i $sample_f5 -r -s $sample_fq -x 'cuda:0,1,2,3' --trim_adapters --bam_out --moves_out
    else
        echo "guppy basecaller not version 6.4, you need that for fast5 file formatting"
    fi
    """
}

process catFq {
    publishDir "./output/${runname}/catFq"
    memory params.mem
    cpus params.t

    input:
    tuple val(runname), path(filepath)

    output:
    tuple val(runname), path('*.fastq')

    script:
    """
    cat ${filepath}/pass/*.fastq > ${runname}.fastq
    """
}

process miniMapping {
    publishDir "./output/${runname}/mapped"
    memory params.mem
    cpus params.t

    input:
    path genome
    tuple val(runname), path(catFq)

    output:	
    tuple val(runname), path('*.sam')
	
    script:
    if (params.mode == "ncrna")
        """
        minimap2 -ax sr -d ${genome.baseName}.mmi ${genome}
        minimap2 -ax sr -uf -t ${params.t} --secondary=no ${genome.baseName}.mmi ${catFq} > ${genome.baseName}.sam
        """        
    else if (params.mode == "cdna")
        """
        minimap2 -ax map-ont -d ${genome.baseName}.mmi ${genome}
        minimap2 -ax map-ont -uf -t ${params.t} --secondary=no ${genome.baseName}.mmi ${catFq} > ${genome.baseName}.sam
        """
}

process sam2bam {
    publishDir "./output/${runname}/bam/"

    input:
    tuple val(runname), path(sam)

    output:
    tuple val(runname), path("*.bam"), path("*.bam.bai")

    script:
    """
    samtools view -Sb ${sam} | samtools sort -o ${sam.baseName}.bam -
    samtools index ${sam.baseName}.bam
    """
}

process nanopolishEventalign {
    publishDir "./output/${runname}/np"
    memory params.mem
    cpus params.t

    input:
    tuple val(runname), val(f5path), path(catFq), path(bamFile), path(baiFile)

    output:
    tuple val(runname), path('eventalign.txt'), path('summary.txt')

    script:
    """
    nanopolish index -d "${f5path}/fast5_pass" ${catFq}
    nanopolish eventalign --reads ${catFq} --bam ${bamFile} --genome ${params.genome} --signal-index --scale-events --summary summary.txt --threads ${params.t} > eventalign.txt
    """
}

process xporeDataprep {
	container 'xpore'
	publishDir "./output/${runname}/${runname}_dataprep", mode: 'copy'
	memory params.mem
	cpus params.t

	input:
	tuple val(runname), path(eventalign), path(summary)

	output: 
	tuple val(runname), val("${params.homeDir}/output/${runname}/${runname}_dataprep")
	file('*')

	shell:
	"""
	xpore dataprep --eventalign $eventalign --n_processes ${params.t} --out_dir "."
	"""
}

process makeYml {
	publishDir "./output/comparisons/${pairnum}/", mode: 'copy'
	memory params.mem
	cpus params.t

	input:
	tuple val(pairnum), val(name1), val(labelA), val(p1), val(name2), val(labelB), val(p2)

	output:
	tuple val(pairnum), file('*.yml'), val(name1), val(name2), val(p1), val(p2)

	shell:
	"""
	${params.ymlScript} -d ${params.ymlDefault} -A ${p1} -B ${p2} > "${name1}_vs_${name2}.yml"
	"""
}

process xporeDiffmod {
	container 'xpore'
	publishDir "./output/comparisons/${pairnum}/", mode: 'copy'
	memory params.mem
	cpus params.t

	input:
	tuple val(pairnum), file(configFile), val(name1), val(name2), path(path1), path(path2)

	output:
	tuple val(pairnum), val("${params.homeDir}/output/comparisons/${pairnum}"), val(name1), val(name2)
	file('*')

	shell:
	"""
	xpore diffmod --config ${configFile} --n_processes ${params.t}
	"""
}

process xporePostprocessing {
	container 'xpore'
	publishDir "./output/comparisons/${pairnum}/", mode: 'copy'
	memory params.mem
	cpus params.t

	input:
	tuple val(pairnum), path(diffmodPath), val(name1), val(name2)

	output:
	tuple val(pairnum), val("${params.homeDir}/output/comparisons/${pairnum}")

	shell:
	"""
	xpore postprocessing --diffmod_dir ${diffmodPath}
	"""
}

process filterXpore {
	conda params.condaenv
	publishDir "./output/comparisons/${pairnum}/", mode: 'copy'
	memory params.mem
	cpus params.t

	input:
	tuple val(pairnum), path(xporePath)

	output:
	tuple val(pairnum), val("${params.homeDir}/output/comparisons/${pairnum}/filtered_xpore.table")
	file('*')

	shell:
	"""
	${params.xfiltScript} --in ${xporePath}/majority_direction_kmer_diffmod.table --out filtered_xpore.table
	"""
}

process xporeReadnames {
	conda params.condaenv
	publishDir "./output/comparisons/${pairnum}/", mode: 'copy'
	memory params.mem
	cpus params.t

	input:
	tuple val(pairnum), val(xporePath), val(Aname), val(A), val(Adataprep), val(Bname), val(B), val(Bdataprep)

	output:
	tuple val(pairnum), file('*')

	shell:
	"""
	python ${params.readnameScript} \
		--xpore_table ${xporePath} \
		--sampleA_dir ${Adataprep} \
		--sampleB_dir ${Bdataprep} \
		--outfile xpore_readnames.csv
	"""
}

process xporeTraindata {
	cache false
	conda params.condaenv
    publishDir "./output/${samplename}/kmer_table"
    memory params.mem
    cpus params.t
    
	input:
	tuple val(runname), file(readnames), val(samplename), val(sample_fq), val(sample_f5)

	output:
	file('*')

	shell:
	"""
	python ${params.tableprepScript} \
		--readname_file ${readnames} \
		--bamdir ${sample_fq} \
		--fdir ${sample_f5} \
		--outfile ${samplename}_table.parquet 
	"""
}
