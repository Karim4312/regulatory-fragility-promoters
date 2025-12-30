import gzip

def extract_tss(gtf_path, output_tsv):
    with open(output_tsv, 'w') as out:
        out.write('gene_id\tgene_name\tchr\tTSS\tstrand\n')

        with gzip.open(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields[2] != 'transcript':
                    continue

                chr_, _, _, start, end, _, strand, _, attributes = fields

                attr = {}
                for a in attributes.strip().split(';'):
                    if ' ' in a:
                        k, v = a.strip().split(' ', 1)
                        attr[k] = v.strip('"')

                gene_id = attr.get('gene_id', '')
                gene_name = attr.get('gene_name', '')
                tss = start if strand == '+' else end

                out.write(f"{gene_id}\t{gene_name}\t{chr_}\t{tss}\t{strand}\n')


if __name__ == "__main__":
    extract_tss(
        "/content/data/gencode.v38.annotation.gtf.gz",
        "/content/data/tss_from_gencode.tsv"
    )
