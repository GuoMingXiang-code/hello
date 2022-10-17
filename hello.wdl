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
   docker: 'ubuntu:impish-20220105'
  }
}

workflow testname {
  meta {
    author: "gmx test"
    email: "gmx18525373056@163.com"
  }
  call hello
}
