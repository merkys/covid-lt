rule S1_antibody_complexes:
    output:
        "S1-antibody-complexes.lst"
    shell:
        "./get-S1-antibody-complexes > {output}"

rule download_pdb:
    output:
        "pdb/{id}.pdb"
    shell:
        "wget https://www.crystallography.net/pdb/{wildcards.id}.pdb -O {output}"
