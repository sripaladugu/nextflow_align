#!/usr/bin/env nextflow

println("Nextflow started")
println(params.fastq)
demuxed = Channel.fromPath(params.fastq)
read_pairs = demuxed
               .map { it -> [ sample(it), it ]}
               .groupTuple( sort: true )
               .map{ sample, reads ->
                 tuple( sample, reads[0], reads[1])}

//
// Step 1. Map the reads using bwa.
// 

genome_ref = file(params.genome_file)
genome_ref_dir = genome_ref.parent
genome_ref_fname = genome_ref.getName()
process mapping {
    publishDir "${params.results_dir}/${sample_id}/", mode: 'copy', overwrite: true

    input:
    set val(sample_id), file(read1), file(read2) from read_pairs
    file genome_ref_dir
    
    output:
    set sample_id, file('alignment.bam') into bam_files

    script:
    """
    echo "${sample_id}"
    bwa mem -M -t 24 ${genome_ref_dir}/${genome_ref_fname} ${read1} ${read2} | samtools view -Sb - > alignment.bam
    """
}

/**
 * @return gets sample_name
 */
def sample(fileproc) {
  def procfile = fileproc.getFileName().toString()
  if(procfile =~ /.+[_\.][Rr][12].+/) {
    value = procfile =~ /(.+)[_\.][Rr][12].+/
  }
  else if(procfile =~ /.+[_\.][12].+/) {
    value = procfile =~ /(.+)[_\.][12].+/
  }
  return value[0][1]
}

/**
 *
 */
 workflow.onComplete {
    println ""
    println "WORKFLOW SUMMARY"
    println "Pipeline started at: $workflow.start"
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Results Directory: ${params.results_dir}"
    println "Error message: $workflow.errorMessage"
    println "Execution status: ${ workflow.success ? 'SUCCESS' : 'failed' }"
}
