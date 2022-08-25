version 1.0

task hello {
  input {
    String name
  }

  command {
    echo 'hello ${name}!'
  }

  output {
    File response = stdout()
  }

  runtime {
        docker: "cr-cn-beijing.volces.com/bio-island/busybox:1.30.0"
        memory: "1 GB"
        cpu: "1"
        bootDiskSizeGb: 10
        disks: "local-disk 40 SSD"
        preemptible: 1
    }
}

workflow testname {
  call hello
}
