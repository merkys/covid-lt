rule get_prodigy_dataset:
    output:
        "dG-datasets/prodigy.csv"
    log:
        "dG-datasets/prodigy.log"
    shell:
        """
        date --rfc-3339=s > {log}
        wget https://bianca.science.uu.nl/prodigy/static/PRODIGY_dataset.csv -O {output} 2>> {log}
        """

rule get_skempi_dataset:
    output:
        "dG-datasets/skempi.csv"
    log:
        "dG-datasets/skempi.log"
    shell:
        """
        date --rfc-3339=s > {log}
        wget https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv -O {output} 2>> {log}
        """

rule merge_dG_datasets:
    input:
        prodigy = "dG-datasets/prodigy.csv",
        skempi = "dG-datasets/skempi.csv"
    output:
        "dG-datasets/merged.tab"
    shell:
        """
        echo PDB chains1 chains2 dG PRODIGY SKEMPI | sed 's/ /\t/g' > {output}
        (
            tail -n +2 {input.prodigy} | cut -d , -f 1,2,4 | awk '{{print $1 ",1,0"}}' | sed 's/.pdb//' | sed 's/:/,/g' | sed 's/,/\t/g'
            tail -n +2 {input.skempi} | cut -d ';' -f 1,9 | sed 's/_/;/g' | sed 's/;/\t/g' \
                | awk '{{print $1 "\t" $2 "\t" $3 "\t" (8.314/4184) * (273.15 + 25.0) * log($4) "\t" 0 "\t" 1}}'
        ) | sort >> {output}
        """
