version development

workflow bitmapper_index {

    input {
        File genome
        String index_name
        String destination
    }
    call index{
                input: genome = genome,
                index_folder = index_name
    }
    call copy {
        input:
            files = [index.out],
            destination = destination
    }

    output {
        Array[File] out = copy.out
    }
}

task index {
   input {
        File genome
        String index_folder
   }

   command {
        /opt/BitMapperBS/bitmapperBS --index ~{genome}  --index_folder ${index_folder}
   }

   output {
        File out = index_folder
   }

  runtime {
    docker: "quay.io/comp-bio-aging/bit_mapper_bs:latest"
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