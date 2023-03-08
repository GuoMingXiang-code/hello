version 1.0

import "./tasks_melt.wdl" as tvc
workflow melt {
    input{
        Array[File] inputBAMs
        Array[File] bamIndexs

        File RefFasta
        File RefIndex
        File RefDict
        File Refamb
        File Refann
        File Refbwt
        File Refpac
        File Refsa

        Array[File] me_refs
        File MELT

    }

    scatter(pair in zip(inputBAMs,bamIndexs)){
        call tvc.call_MELT_step1{
            input:inputBAM=pair.left,
                bamIndex=pair.right,
                me_refs=me_refs,
                RefFasta=RefFasta,
                RefIndex=RefIndex,
                RefDict=RefDict,
                Refamb=Refamb,
                Refann=Refann,
                Refbwt=Refbwt,
                Refpac=Refpac,
                Refsa=Refsa,
                MELT=MELT
        }
    }

    call tvc.call_MELT_step2{
        input: step1_output=flatten(call_MELT_step1.output_step1),
                me_refs=me_refs,
                RefFasta=RefFasta,
                RefIndex=RefIndex,
                RefDict=RefDict,
                Refamb=Refamb,
                Refann=Refann,
                Refbwt=Refbwt,
                Refpac=Refpac,
                Refsa=Refsa,
                MELT=MELT

    }

    scatter(pair in zip(inputBAMs,bamIndexs)){
        call tvc.call_MELT_step3{
            input: step2_output=call_MELT_step2.output_step2,
                inputBAM=pair.left,
                bamIndex=pair.right,
                me_refs=me_refs,
                RefFasta=RefFasta,
                RefIndex=RefIndex,
                RefDict=RefDict,
                Refamb=Refamb,
                Refann=Refann,
                Refbwt=Refbwt,
                Refpac=Refpac,
                Refsa=Refsa,
                MELT=MELT
        }
    }

    call tvc.call_MELT_step4{
        input: step2_output=call_MELT_step2.output_step2,
                step3_output=flatten(call_MELT_step3.output_step3),
                me_refs=me_refs,
                RefFasta=RefFasta,
                RefIndex=RefIndex,
                RefDict=RefDict,
                Refamb=Refamb,
                Refann=Refann,
                Refbwt=Refbwt,
                Refpac=Refpac,
                Refsa=Refsa,
                MELT=MELT
    }
}
