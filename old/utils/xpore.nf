#!/usr/bin/env nextflow

/////////////////////////////////////////////////////////////////
///////////////////////// setup /////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////// setup the sample data ////////////////////////
//// MUST BE nested array where the 1st element is codename and second is path to the parent folder 

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

fh = new File(params.inputFile)
def csv_content = fh.getText('utf-8')
 
def data_iterator = parseCsv(csv_content, separator: '\t', readFirstLine: false)
// println data_iterator.getClass()  // class com.xlson.groovycsv.CsvIterator
 
def sampleArray = []
def signalArray = []
def sampleArrayPairs = []
i = 0
//Sample-Codename FqPath  f5Path  ComparisonPair  ComparisonOrder
for (line in data_iterator) {
    sampleArray[i] = [line[0],line[1]]
    signalArray[i] = [line[0],line[2]]
    sampleArrayPairs[i] = [line[0],line[3],line[4]]
    i += 1;
}

runDataChannel = Channel.from(sampleArray).unique()
signalChannel = Channel.from(signalArray).unique()
	
////////////////////////////////////////////////////////////
//////////////////////// dataprep //////////////////////////
////////////////////////////////////////////////////////////

// once for N samples (they use the same index)
if (params.mode == "cdna") {
  process genomeIndexCdna {
    publishDir "./output/genomeIndexed"
    memory params.mem
    cpus params.t

    output:
    file('*.mmi') into minimapIndex

    script:
    """
    minimap2 -ax map-ont -d genome_index.mmi ${params.genome}
    """
  }
}
if (params.mode == "ncrna") {
  process genomeIndexNcrna {
    publishDir "./output/genomeIndexed"
    memory params.mem
    cpus params.t

    output:
    file('*.mmi') into minimapIndex

    script:
    """
    minimap2 -ax sr -d genome_index.mmi ${params.genome}
    """
  }
}

// once per sample
process catfq {
  publishDir "./output/$runname/catFq"

  memory params.mem
  cpus params.t

  input:
  tuple runname, filepath from runDataChannel

  output:
  tuple runname, '*.fastq' into joinFq
  tuple runname, '*.fastq' into joinFqMirror

  script:
  """
  cat ${filepath}/pass/*.fastq > ${runname}.fastq
  """
}

minimapChan = joinFq.combine(minimapIndex)

if (params.mode == "cdna") {
  process miniMappingCdna {
    publishDir "./output/$runname/mapped"
    memory params.mem
    cpus params.t

    input:
    tuple runname, catFq, indexFile from minimapChan

    output:
    tuple runname, '*.sam' into samFiles

    script:
    """
    minimap2 -ax map-ont -uf -t ${params.t} --secondary=no ${indexFile} ${catFq} > ${indexFile.baseName}.sam
    """
  }
}
if (params.mode == "ncrna") {
  process miniMappingNcrna {
    publishDir "./output/$runname/mapped"
    memory params.mem
    cpus params.t

    input:
    tuple runname, catFq, indexFile from minimapChan

    output:
    tuple runname, '*.sam' into samFiles

    script:
    """
    minimap2 -ax sr -uf -t ${params.t} --secondary=no ${indexFile} ${catFq} > ${indexFile.baseName}.sam
    """
  }
}

process sam2bam {
  publishDir "./output/$runname/bam"
  memory params.mem
  cpus params.t

  input:
  tuple runname, samFile from samFiles

  output:
  tuple runname, '*.bam' into bamFiles

  script:
  """
  samtools view -Sb ${samFile} | samtools sort -o ${samFile.baseName}.bam -
  samtools index ${samFile.baseName}.bam
  """
}

///////// combine using shared keys (runname) ////////////////////
npInputChan = signalChannel.join(joinFqMirror).join(bamFiles)

process nanopolish {
  publishDir "./output/$runname/np"
  memory params.mem
  cpus params.t

  input:
  tuple runname, filepath, catFq, bamFile from npInputChan

  output:
  tuple runname, 'eventalign.txt' into eventAlign
  tuple runname, 'summary.txt' into summary

  script:
  """
  nanopolish index -d "${filepath}/fast5_pass" ${catFq}
  nanopolish eventalign --reads ${catFq} --bam ${bamFile} --genome ${params.genome} --signal-index --scale-events --summary summary.txt --threads ${params.t} > eventalign.txt
  """
}

// NOTE: outputs look weird because we need to pass one path which gets read as a path, and one which gets read as a value
process dataprep {
  container 'xpore'
  publishDir "./output/$runname/${runname}_dataprep", mode: 'copy'
  memory params.mem
  cpus params.t

  input:
  tuple val(runname), path(eventalign) from eventAlign

  output: 
  file('*') into xporeOut
  tuple val(runname), val("${params.homeDir}/output/${runname}/${runname}_dataprep") into dataprep

  shell:
  """
  xpore dataprep --eventalign $eventalign --n_processes ${params.t} --out_dir "."
  """
}

//////////////////////////////////////////////////////
/////////////////////// DIFFMOD //////////////////////
//////////////////////////////////////////////////////

// preprocess for diffmod, create YML files for comparisons
// NOTE:: we are doubling up on the paths, here's why
// when I use path(VAR) as INPUT, it enables processes to get the local, ie work dir, links to the files
// skipping this results in processes which lack the required files
// HOWEVER, the diffmod process requires access to the dataprep directories, but doesn't take them directly as input
// so what I do is double-up the paths in the OUTPUT of dataprep, so I can use ONE copy as path(v), and another as val(v)
// then that val(v) can go thru makeYml, and become path(v) in the diffmod
// slightly annoying

if (params.pairwise) {

  // dataprep: SAMPLE, PATH   
  dataprep.into { dp1; dp2 }
  basedata = Channel.from(sampleArrayPairs)
  sampleA = Channel.create()
  sampleB = Channel.create()
  
  // sampleA: SAMPLE, NUMBER, A   
  basedata.choice( sampleA, sampleB ) { a -> a[2] =~ 'A' ? 0 : 1 }
  // As: SAMPLE, NUMBER, A, PATH
  As = sampleA.combine(dp1, by: 0)
  Bs = sampleB.combine(dp2, by: 0)
  // pairChannel: NUMBER, SAMPLE, AB, PATH, SAMPLE, AB, PATH
  pairChannel = As.combine(Bs, by: 1)

  process makeYml {
    publishDir "./output/comparisons/$pairnum/", mode: 'copy'
    memory params.mem
    cpus params.t

    input:
    tuple val(pairnum), val(name1), val(labelA), val(p1), val(name2), val(labelB), val(p2) from pairChannel

    output:
    tuple val(pairnum), file('*.yml'), val(name1), val(name2), val(p1), val(p2) into diffmodInput

    shell:
    """
    ${params.ymlScript} -d ${params.ymlDefault} -A ${p1} -B ${p2} > "${name1}_vs_${name2}.yml"
    """
  }

  process diffmod {
    container 'xpore'
    publishDir "./output/comparisons/$pairnum/", mode: 'copy'
    memory params.mem
    cpus params.t

    input:
    tuple val(pairnum), file(configFile), val(name1), val(name2), path(path1), path(path2) from diffmodInput

    output:
    file('*') into diffmodOut
    tuple val(pairnum), val("${params.homeDir}/output/comparisons/$pairnum"), val(name1), val(name2) into ppChan

    shell:
    """
    xpore diffmod --config ${configFile} --n_processes ${params.t}
    """
  }

  ppChan.into { xporePP; addchrPP }

  process postprocessing {
    container 'xpore'
    publishDir "./output/comparisons/$pairnum/", mode: 'copy'
    memory params.mem
    cpus params.t

    input:
    tuple val(pairnum), path(diffmodPath), val(name1), val(name2) from xporePP

    shell:
    """
    xpore postprocessing --diffmod_dir ${diffmodPath}
    """
  }
}
