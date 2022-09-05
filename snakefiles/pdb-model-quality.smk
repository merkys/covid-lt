rule voromqa:
    input:
        "{path}/{pdbid}.pdb"
    output:
        "{path}/{pdbid}.voromqa"
    singularity:
        "container.sif"
    shell:
        "voronota-voromqa -i {input} | cut -d ' ' -f 2- > {output} || echo WARNING: {output}: rule failed >&2"

def pristine_pdbs_voromqa(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(pdb_input_dir + "{pdbid}.voromqa", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def fixed_pdbs_voromqa(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/fixed/{pdbid}.voromqa", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def optimized_pdbs_voromqa(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/optimized/{pdbid}.voromqa", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule voromqa_all:
    input:
        pristine_pdbs_voromqa = pristine_pdbs_voromqa,
        fixed_pdbs_voromqa = fixed_pdbs_voromqa,
        optimized_pdbs_voromqa = optimized_pdbs_voromqa
    output:
        output_dir + "voromqa.tab"
    shell:
        """
        for VOROMQA in {input.pristine_pdbs_voromqa}
        do
            PRISTINE=$VOROMQA
            FIXED={output_dir}pdb/fixed/$(basename $VOROMQA)
            OPTIMIZED={output_dir}pdb/optimized/$(basename $VOROMQA)
            if [ -s $PRISTINE -a -s $FIXED ]
            then
                (
                    echo -en $(basename $VOROMQA .voromqa)"\t"
                    paste $PRISTINE $FIXED $OPTIMIZED
                ) | xargs echo | sed 's/ /\t/g'
            fi
        done > {output}
        """

# FIXME: bin/qmean depends on qmean Python module, not yet in Debian.
rule qmean:
    input:
        "{path}/{pdbid}.pdb"
    output:
        "{path}/{pdbid}.qmean"
    log:
        "{path}/{pdbid}.qmean.log"
    singularity:
        "container.sif"
    shell:
        """
        if ! bin/qmean {input} > {output} 2> {log}
        then
            echo WARNING: {output}: rule failed >&2
            cat {log} >&2
        fi
        """

def pristine_pdbs_qmean(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(pdb_input_dir + "{pdbid}.qmean", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def fixed_pdbs_qmean(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/fixed/{pdbid}.qmean", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

def optimized_pdbs_qmean(wildcards):
    from glob import glob
    checkpoint_output = checkpoints.download_pdb_all.get(**wildcards).output[0]
    return expand(output_dir + "pdb/optimized/{pdbid}.qmean", pdbid=glob_wildcards(checkpoint_output + '/{pdbid}.pdb').pdbid)

rule qmean_all:
    input:
        pristine_pdbs_qmean = pristine_pdbs_qmean,
        fixed_pdbs_qmean = fixed_pdbs_qmean,
        optimized_pdbs_qmean = optimized_pdbs_qmean
    output:
        output_dir + "qmean.tab"
    shell:
        """
        for QMEAN in {input.pristine_pdbs_qmean}
        do
            PRISTINE=$QMEAN
            FIXED={output_dir}pdb/fixed/$(basename $QMEAN)
            OPTIMIZED={output_dir}pdb/optimized/$(basename $QMEAN)
            if [ -s $PRISTINE -a -s $FIXED ]
            then
                (
                    echo -en $(basename $QMEAN .QMEAN)"\t"
                    paste $PRISTINE $FIXED $OPTIMIZED
                ) | xargs echo | sed 's/ /\t/g'
            fi
        done > {output}
        """
