#!/usr/bin/env python3
"""
从基因名提取所有转录本的uORF序列并输出FASTA格式
支持经典和非经典起始密码子
"""

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import argparse

# 设置你的email（NCBI要求）
Entrez.email = "your.email@example.com"


def get_transcripts_from_gene(gene_name, species="human"):
    """
    从基因名获取所有转录本序列

    参数:
    gene_name: 基因名（如TP53, BRCA1）
    species: 物种（human, mouse等）

    返回:
    转录本列表，每个包含序列和CDS位置信息
    """
    print(f"正在搜索基因: {gene_name} ({species})...", file=sys.stderr)

    # 构建搜索词
    if species.lower() == "human":
        search_term = f"{gene_name}[Gene Name] AND Homo sapiens[Organism] AND RefSeq[Filter]"
    elif species.lower() == "mouse":
        search_term = f"{gene_name}[Gene Name] AND Mus musculus[Organism] AND RefSeq[Filter]"
    else:
        search_term = f"{gene_name}[Gene Name] AND {species}[Organism] AND RefSeq[Filter]"

    # 搜索基因对应的核酸序列
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=100)
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        print(f"找到 {len(id_list)} 个转录本记录", file=sys.stderr)

        if not id_list:
            print(f"警告: 未找到基因 {gene_name} 的转录本", file=sys.stderr)
            return []

        # 获取序列详细信息
        handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
        records = SeqIO.parse(handle, "genbank")

        transcripts = []
        for rec in records:
            # 提取CDS位置
            cds_start = None
            for feature in rec.features:
                if feature.type == "CDS":
                    # 获取CDS起始位置（第一个外显子的起始）
                    cds_start = int(feature.location.start)
                    break

            if cds_start is not None:
                transcripts.append({
                    'id': rec.id,
                    'description': rec.description,
                    'sequence': str(rec.seq),
                    'cds_start': cds_start
                })
                print(f"  - {rec.id}: CDS起始于位置 {cds_start}", file=sys.stderr)

        handle.close()
        return transcripts

    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        return []


def find_all_uorfs(sequence, cds_start, include_noncanonical=True, min_length=9):
    """
    搜索所有uORF（包括非经典起始密码子）
    """
    sequence = sequence.upper()

    # 定义起始和终止密码子
    if include_noncanonical:
        start_codons = ['ATG', 'CTG', 'GTG', 'TTG', 'ACG', 'ATT', 'ATA']
    else:
        start_codons = ['ATG']

    stop_codons = ['TAA', 'TAG', 'TGA']
    uorfs = []

    # 在5' UTR中搜索（从0到cds_start之前）
    for start_pos in range(cds_start - 2):
        codon = sequence[start_pos:start_pos + 3]

        if len(codon) < 3:
            continue

        if codon in start_codons:
            # 确定阅读框
            frame = start_pos % 3

            # 在该阅读框中延伸寻找终止密码子
            found_stop = False
            for stop_pos in range(start_pos + 3, len(sequence), 3):
                # 保持在同一阅读框
                if stop_pos % 3 != frame:
                    continue

                stop_codon = sequence[stop_pos:stop_pos + 3]

                if len(stop_codon) < 3:
                    break

                if stop_codon in stop_codons:
                    uorf_length = stop_pos + 3 - start_pos

                    # 过滤最小长度
                    if uorf_length >= min_length:
                        # 判断类型
                        if stop_pos + 3 <= cds_start:
                            uorf_type = 'independent'
                        else:
                            uorf_type = 'overlapping'

                        uorfs.append({
                            'start': start_pos,
                            'end': stop_pos + 3,
                            'start_codon': codon,
                            'stop_codon': stop_codon,
                            'length': uorf_length,
                            'frame': frame,
                            'type': uorf_type,
                            'sequence': sequence[start_pos:stop_pos + 3]
                        })

                    found_stop = True
                    break

            # 如果没有终止密码子，延伸到CDS起始或序列末端
            if not found_stop:
                # 计算到CDS或序列末端的距离
                end_pos = min(cds_start, len(sequence))
                # 调整到阅读框对齐的位置
                while end_pos > start_pos + 3 and (end_pos - start_pos) % 3 != 0:
                    end_pos -= 1

                if end_pos > start_pos + 3:
                    uorf_length = end_pos - start_pos
                    if uorf_length >= min_length:
                        uorfs.append({
                            'start': start_pos,
                            'end': end_pos,
                            'start_codon': codon,
                            'stop_codon': 'NONE',
                            'length': uorf_length,
                            'frame': frame,
                            'type': 'no_stop',
                            'sequence': sequence[start_pos:end_pos]
                        })

    return uorfs


def uorfs_to_fasta(gene_name, transcripts, include_noncanonical=True, min_length=9):
    """
    将所有转录本的uORF转换为FASTA格式
    """
    fasta_records = []
    bed_records = []
    total_uorfs = 0

    for transcript in transcripts:
        transcript_id = transcript['id']
        sequence = transcript['sequence']
        cds_start = transcript['cds_start']

        # 查找uORF
        uorfs = find_all_uorfs(sequence, cds_start, include_noncanonical, min_length)

        print(f"\n转录本 {transcript_id}: 找到 {len(uorfs)} 个uORF", file=sys.stderr)

        for idx, uorf in enumerate(uorfs, 1):
            # 构建FASTA header
            header = (
                f"{gene_name}|{transcript_id}|uORF{idx}|"
                f"{uorf['start_codon']}|"
                f"pos:{uorf['start']}-{uorf['end']}|"
                f"len:{uorf['length']}bp|"
                f"frame:{uorf['frame']}|"
                f"type:{uorf['type']}"
            )

            # 创建SeqRecord
            seq_record = SeqRecord(
                Seq(uorf['sequence']),
                id=header,
                description=""
            )

            fasta_records.append(seq_record)

            # 创建BED记录
            bed_name = f"{gene_name}_uORF{total_uorfs + 1}"
            bed_record = {
                'chrom': transcript_id,
                'start': uorf['start'],
                'end': uorf['end'],
                'name': bed_name,
                'score': uorf['length'],  # 使用长度作为score
                'strand': '+',  # 转录本坐标系始终是正向
                'start_codon': uorf['start_codon'],
                'stop_codon': uorf['stop_codon'],
                'type': uorf['type'],
                'frame': uorf['frame']
            }
            bed_records.append(bed_record)

            total_uorfs += 1

    print(f"\n总计: {total_uorfs} 个uORF", file=sys.stderr)
    return fasta_records, bed_records


def write_bed_file(bed_records, output_file):
    """
    写入BED格式文件（BED12格式，包含额外信息）
    """
    with open(output_file, 'w') as f:
        # 写入header（可选，便于理解）
        f.write("# chrom\tstart\tend\tname\tscore\tstrand\tstart_codon\tstop_codon\ttype\tframe\n")

        for record in bed_records:
            line = (
                f"{record['chrom']}\t"
                f"{record['start']}\t"
                f"{record['end']}\t"
                f"{record['name']}\t"
                f"{record['score']}\t"
                f"{record['strand']}\t"
                f"{record['start_codon']}\t"
                f"{record['stop_codon']}\t"
                f"{record['type']}\t"
                f"{record['frame']}\n"
            )
            f.write(line)


def main():
    parser = argparse.ArgumentParser(
        description='从基因名提取所有转录本的uORF序列',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s TP53 -o tp53_uorfs.fasta
  %(prog)s BRCA1 --species human --canonical-only
  %(prog)s Actb --species mouse --min-length 15
        """
    )

    parser.add_argument('gene', help='基因名（如TP53, BRCA1）')
    parser.add_argument('-o', '--output', default='uorfs.fasta',
                        help='输出FASTA文件名（默认: uorfs.fasta）')
    parser.add_argument('-b', '--bed', default=None,
                        help='输出BED文件名（默认: 与FASTA同名但扩展名为.bed）')
    parser.add_argument('-s', '--species', default='human',
                        help='物种（默认: human）')
    parser.add_argument('-c', '--canonical-only', action='store_true',
                        help='仅搜索ATG起始的uORF（默认包含非经典起始）')
    parser.add_argument('-m', '--min-length', type=int, default=9,
                        help='最小uORF长度（bp，默认: 9）')
    parser.add_argument('-e', '--email', default='your.email@example.com',
                        help='NCBI Entrez email（必填）')

    args = parser.parse_args()

    # 设置email
    Entrez.email = args.email

    # 获取转录本
    transcripts = get_transcripts_from_gene(args.gene, args.species)

    if not transcripts:
        print("未找到转录本，退出。", file=sys.stderr)
        sys.exit(1)

    # 提取uORF并转换为FASTA
    include_noncanonical = not args.canonical_only
    fasta_records, bed_records = uorfs_to_fasta(
        args.gene,
        transcripts,
        include_noncanonical,
        args.min_length
    )

    # 写入FASTA文件
    if fasta_records:
        SeqIO.write(fasta_records, args.output, "fasta")
        print(f"\n结果已保存到: {args.output}", file=sys.stderr)

        # 写入BED文件
        bed_output = args.bed if args.bed else args.output.rsplit('.', 1)[0] + '.bed'
        write_bed_file(bed_records, bed_output)
        print(f"BED文件已保存到: {bed_output}", file=sys.stderr)
    else:
        print("\n未找到任何uORF", file=sys.stderr)


if __name__ == "__main__":
    main()
