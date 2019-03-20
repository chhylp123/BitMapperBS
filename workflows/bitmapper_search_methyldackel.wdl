version development

workflow bs_map {
    input {
        File genome_index
        Array[File] reads
        String filename
        String destination
        Boolean is_paired = true
        Int map_threads = 8
    }

    call bitmapper {
        input:
            index_folder = genome_index,
            reads = reads,
            is_paired = is_paired,
            filename = filename,
            threads = map_threads
    }

    call copy as copy_bit {
        input:
            files = [bitmapper.out, bitmapper.stats],
            destination = destination
    }

    call picard_readgroups_sort {
        input: bam = bitmapper.out,
                    filename = filename
    }

    call copy as copy_sorted {
        input:
        files = [picard_readgroups_sort.out],
        destination = destination
    }


}

task bitmapper {
   input {
        File index_folder
        Array[File] reads
        Boolean is_paired
        String filename
        Int threads
   }
   command {
        /opt/BitMapperBS/bitmapperBS --search ~{index_folder} ~{if(is_paired) then " --seq1 " + reads[0] + " --seq2 "+ reads[1] + " --sensitive --pe" else " --seq1 " + reads[0]} -t ~{threads} --mapstats --bam -o ~{filename}.bam
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:latest"
  }

  output {
    File out = "~{filename}.bam"
    File stats = "~{filename}.bam.mapstats"
  }
}

task picard_readgroups_sort{
    input {
        File bam
        String filename
    }
    command {
        picard AddOrReplaceReadGroups \
        I=~{bam} \
        O=~{filename}_sorted.bam \
        RGID=4 \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=20 \
        SORT_ORDER=coordinate
    }

    runtime {
        docker: "biocontainers/picard:v2.3.0_cv3"
    }

    output {
        File out = "~{filename}_sorted.bam"
    }

}

task methyldackel {
    input {
        String name
        File bam
        File genome
        Int threads = 4
    }

    command {
        MethylDackel extract --CHH --CHG --counts -@ ~{threads} ~{genome} ~{bam}
    }


    runtime {
        docker: "quay.io/biocontainers/methyldackel@sha256:d434c3e320a40648a3c74e268f410c57649ab208fcde4da93677243b22900c55" #0.3.0--h84994c4_3
    }

    output {
       File chg = name + "_CHG.counts.bedGraph"
       File chh = name + "_CHH.counts.bedGraph"
       File cpg = name + "_CpG.counts.bedGraph"
    }
}

task copy {
    input {
        Array[File] files
        String destination
    }

    command {
        mkdir -p ~{destination}
        cp -L -R -u ~{sep=' ' files} ~{destination}
    }

    output {
        Array[File] out = files
    }
}