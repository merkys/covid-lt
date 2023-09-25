rule R_images:
    input:
        "doc/images/{base}.R"
    output:
        "doc/images/{base}.svg"
    shell:
        "{input} > {output}"
