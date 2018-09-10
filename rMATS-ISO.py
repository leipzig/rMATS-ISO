#!/usr/bin/env python3

import sys
import re
import argparse
import utils
import os


class iso_para:
    def __init__(self, args):
        # input file
        ## genome and annotation
        self.in_genome = os.path.abspath(args.in_genome) if args.in_genome else ''
        self.star_idx = os.path.abspath(args.STAR_idx) if args.STAR_idx else self.in_genome + '.STAR'
        self.mmap2_idx = os.path.abspath(args.minimap2_idx) if args.minimap2_idx else self.in_genome + '.smmi'
        self.in_gtf = os.path.abspath(args.in_gtf) if args.in_gtf else ''
        self.rm_gtf = os.path.abspath(args.rm_gtf) if args.rm_gtf else ''
        ## BAM and FASTA
        self.in_bam_list = args.in_bam
        self.in_bam_str = ''
        self.in_long_list = args.in_long
        self.in_short_list = args.in_short

        # output director
        self.out_dir = os.path.abspath(args.out_dir)
        self.updated_gtf = self.out_dir + '/lr2rmats/output/updated.gtf'

        # parameters
        self.nthread = args.thread
        ## lr2rmats
        self.long_aln_cov = args.aln_cov
        self.long_iden_frac = args.iden_frac
        self.long_sec_rat = args.sec_rat
        self.long_sup_cnt = args.sup_cnt
        self.long_split_trans = args.split_trans

        ## ISO_module
        self.module_novel_sj = '-n' if args.novel_sj else ''
        self.module_use_un_pair = '-u' if args.use_un_pair else ''
        self.module_use_multi = '-m' if args.use_multi else ''
        self.module_intron_len = args.intron_len
        self.module_exon_thres = args.exon_thres
        self.module_junc_cnt = args.junc_cnt
        self.module_novel_jcnt = args.novel_jcnt
        self.module_mdl_exon_thres = args.mdl_exon_thres
        self.module_iso_cnt_thres = args.iso_cnt_thres
        self.module_recursive = '-r' if args.recursive else ''
        self.module_iso_exon = ''
        self.module_iso_matrix1 = ''
        self.module_iso_matrix2 = ''

        ## rmats-em
        self.em_seed = args.random_seed
        self.em_n_iter = args.iter_num
        self.em_draw_samp = args.draw_samp
        self.em_pair_simu = args.pair_simu
        self.em_epsilon = args.epsilon
        self.em_p_thres = args.p_thres
        self.EM_out = ''

        ## iso-classify && iso-plot
        self.module_type_summary = ''
        self.module_type = ''
        self.module_coor = ''
        self.module_gene = ''


def parse_argv(parser):
    parser.add_argument('--in-gtf', type=str, required=True, help='Gene annotation in GTF format.')
    read_argv = parser.add_mutually_exclusive_group(required=True)
    read_argv.add_argument('--in-bam', type=str, default='',
                           help='Short-read BAM input file list. Use \'--in-short\' if no BAM file is available. Refer to \'short_read_bam_input.list.example\'')
    read_argv.add_argument('--in-short', type=str, default='',
                           help='Short-read input file list. Refer to \'short_read_fa_input.list.example\'')
    parser.add_argument('-o', '--out-dir', type=str, required=True, help='Output directory.')
    parser.add_argument('-t', '--thread', type=int, default=4, help='Number of threads to use.')

    input_argv = parser.add_argument_group('Optional input files')
    input_argv.add_argument('--in-genome', type=str, default='',
                            help='Genome sequence file in FASTA format. When \'--in-short\' or \'--in-long\' is used, genome is used to do alignment.')
    input_argv.add_argument('--STAR-idx', type=str, default='',
                            help='STAR index directory for genome sequence. Index will be built if not set.')
    input_argv.add_argument('--minimap2-idx', type=str, default='',
                            help='Minimap2 index directory for genome sequence. Index will be built if not set.')
    input_argv.add_argument('--in-long', type=str,
                            help='Long-read input file list. Updated gene annotation will be used when long-read data is provided. Refer to \'long_read_fa_input.list.example\'')
    input_argv.add_argument('--rm-gtf', type=str, default='',
                            help='Transcript file in GTF format that need to be removed from long-read data, e.g., rRNA.')

    long_argv = parser.add_argument_group('Generation of updated annotation from long-read data')
    long_argv.add_argument('--aln-cov', type=float, default=0.67,
                           help='Minimum fraction of long-read\'s aligned bases.')
    long_argv.add_argument('--iden-frac', type=float, default=0.75,
                           help='Minimum fraction of long-read\'s identically aligned bases.')
    long_argv.add_argument('--sec-rat', type=float, default=0.98,
                           help='Maximum ratio of second best and best score to retain the best alignment.')
    long_argv.add_argument('--sup-cnt', type=int, default=1,
                           help='Minimum supporting count of novel-junction spanning short-read to make novel-junction reliable.')
    long_argv.add_argument('--split-trans', default=False, action='store_true',
                           help='Split long-read derived transcript on unreliable junctions.')

    module_argv = parser.add_argument_group('Splicing module detection and quantification')
    module_argv.add_argument('--novel-sj', action='store_true', default=False,
                             help='Allow novel splice junction in the ASM.\nNovel splice junction here means novel combination of known splice sites.')
    module_argv.add_argument('--use-un-pair', action='store_true', default=False,
                             help='Use both proper and improper paired mapped reads.')
    module_argv.add_argument('--use-multi', action='store_true', default=False,
                             help='Use both uniq- and multi-mapped reads.')
    module_argv.add_argument('--fil-mul-gene', action='store_true', default=False,
                             help='Filter out module that overlapps with multi-gene region.')
    module_argv.add_argument('--intron-len', type=int, default=3,
                             help='Minimum intron length for junction read.')
    module_argv.add_argument('--exon-thres', type=int, default=10,
                             help='Maximum number of exon for ASM to enumerate all possible isoforms.\nFor larger ASM, heuristic strategy will be applied.')
    module_argv.add_argument('--junc-cnt', type=int, default=1,
                             help='Minimum allowed number of read count for known-junction.')
    module_argv.add_argument('--novel-jcnt', type=int, default=3,
                             help='Minimum allowed number of read count for novel-junction.')
    module_argv.add_argument('--mdl-exon-thres', type=int, default=50,
                             help='Maximum allowed number of exons in one ASM.')
    module_argv.add_argument('--iso-cnt-thres', type=int, default=20,
                             help='Maximum allowed number of isoforms in one ASM.')
    module_argv.add_argument('--recursive', action='store_true', default=False,
                             help='Recursively output alternative splicing module which is fully embedded in a larger module.')

    em_argv = parser.add_argument_group('Statistical test of splicing difference')
    em_argv.add_argument('--random-seed', type=int, default=17495,
                         help='')
    em_argv.add_argument('--draw-samp', type=int, default=500,
                         help='Number of samples to draw for importance sampling.')
    em_argv.add_argument('--iter-num', type=int, default=100,
                         help='Maximum number of EM iterations to perform.')
    em_argv.add_argument('--pair-simu', type=int, default=1e+05,
                         help='Number of simulations to perform to compute pairwise isoform p-values.')
    em_argv.add_argument('--epsilon', type=float, default=0.01,
                         help='Error tolerance for terminating the EM algorithm.')
    em_argv.add_argument('--p-thres', type=float, default=1.0,
                         help='Below which module-wise p-value shall individual isoform differences be tested.')

    # type_argv = parser.add_argument_group('Classification of types of alternative splicing:\n')

    # visual_argv = parser.add_argument_group('Visualization using sashimi plots.:\n')


def lr2rmats(para, cur_dir):
    utils.exec_cmd(sys.stderr, 'lr2rmats', 'mkdir {}/lr2rmats 2> /dev/null'.format(para.out_dir))
    star_idx = '--STAR-idx {} '.format(para.star_idx) if para.star_idx else ''
    mmap2_idx = '--minimap2-idx {} '.format(para.mmap2_idx) if para.mmap2_idx else ''
    rm_gtf = '--rm-gtf {} '.format(para.rm_gtf) if para.rm_gtf else ''
    split_trans = '--split-trans ' if para.long_split_trans else ''
    op = star_idx + mmap2_idx + rm_gtf + split_trans
    cmd = 'python {}/lr2rmats/run_snakemake.py --cores {} --genome {} --gtf {} --long-read {} --short-read {} --aln-cov {} --iden-frac {} --sec-rat {} --sup-cnt {} --out-dir {}/lr2rmats {}'.format(
        cur_dir, para.nthread, para.in_genome, para.in_gtf, para.in_long_list, para.in_short_list, para.long_aln_cov,
        para.long_iden_frac, para.long_sec_rat, para.long_sup_cnt, para.out_dir, op)
    utils.exec_cmd(sys.stderr, 'lr2rmats', cmd)
    return


def short_map(para):
    utils.exec_cmd(sys.stderr, 'STAR_mapping', 'mkdir {}/STAR_mapping 2> /dev/null'.format(para.out_dir))
    para.in_bam_str = ''
    if not os.path.exists(para.star_idx):
        idx_log = '{}/STAR_mapping/idx.log'.format(para.out_dir)
        cmd = 'STAR --runMode genomeGenerate --runThreadN {} --genomeFastaFiles {} --genomeDir {} --outFileNamePrefix {} >> {}'.format(
            para.nthread, para.in_genome, para.star_idx, idx_log, idx_log)
        utils.exec_cmd(sys.stderr, 'STAR_indexing', cmd)
    with open(para.in_short_list) as in_list:
        line = in_list.readline()
        line = line[:line.index('#')] if '#' in line else line
        n_samp_short = int(line.split()[0])
        samp_i = 1
        for i in range(0, n_samp_short):
            in_bam = []
            line = in_list.readline()
            line = line[:line.index('#')] if '#' in line else line
            n_rep = int(line.split()[0])
            for j in range(0, n_rep):
                line = in_list.readline()
                line = line[:line.index('#')] if '#' in line else line
                in_read = line[:-1]
                log = '{}/STAR_mapping/samp{}.STAR.log'.format(para.out_dir, samp_i)
                cmd = "STAR --runThreadN {}  --genomeDir {} --readFilesIn {} --outFileNamePrefix {}/STAR_mapping/samp{}.STAR --outSAMtype BAM SortedByCoordinate " \
                      "--outFilterType BySJout   --outFilterMultimapNmax 20 --outFilterMismatchNmax 999   --alignIntronMin 25   " \
                      "--alignIntronMax 1000000   --alignMatesGapMax 1000000 --alignSJoverhangMin 8   --alignSJDBoverhangMin 5  " \
                      " --sjdbGTFfile {}  --sjdbOverhang 100 > {} 2 >& 1; samtools index {}/STAR_mapping/samp{}.STARAligned.sortedByCoord.out.bam".format(
                    para.nthread, para.star_idx, in_read, para.out_dir, samp_i, para.in_gtf, log, para.out_dir, samp_i)

                utils.exec_cmd(sys.stderr, 'STAR_mapping', cmd)
                in_bam.append('{}/STAR_mapping/samp{}.STARAligned.sortedByCoord.out.bam'.format(para.out_dir, samp_i))
                samp_i += 1
            para.in_bam_str += ','.join(in_bam) if para.in_bam_str == '' else ':' + ','.join(in_bam)
    para.in_bam_list = ''
    return


def ISO_module(para, cur_dir):
    if para.in_bam_list:  # get bam str from list file
        para.in_bam_str = ''
        with open(para.in_bam_list) as in_list:
            line = in_list.readline()
            line = line[:line.index('#')] if '#' in line else line
            n_samp_short = int(line.split()[0])
            for i in range(0, n_samp_short):
                in_bam = []
                line = in_list.readline()
                line = line[:line.index('#')] if '#' in line else line
                n_rep = int(line.split()[0])
                for j in range(0, n_rep):
                    line = in_list.readline()
                    line = line[:line.index('#')] if '#' in line else line
                    in_bam.append(os.path.abspath(line.split()[0]))
                para.in_bam_str += ','.join(in_bam) if para.in_bam_str == '' else ':' + ','.join(in_bam)

    utils.exec_cmd(sys.stderr, 'IsoModule', 'mkdir {}/ISO_module 2> /dev/null'.format(para.out_dir))
    op = ' '.join([para.module_novel_sj, para.module_use_un_pair, para.module_use_multi, para.module_recursive])
    cmd = 'IsoModule {} {} {} -i {} -T {} -c {} -v {} -e {} -C {} -o {}/ISO_module'.format(             para.in_gtf,
                                                                                                        para.in_bam_str,
                                                                                                        op,
                                                                                                        para.module_intron_len,
                                                                                                        para.module_exon_thres,
                                                                                                        para.module_junc_cnt,
                                                                                                        para.module_novel_jcnt,
                                                                                                        para.module_mdl_exon_thres,
                                                                                                        para.module_iso_cnt_thres,
                                                                                                        para.out_dir)
    utils.exec_cmd(sys.stderr, 'IsoModule', cmd)
    [samp1, samp2] = para.in_bam_str.split(':')
    iso_matrix1 = []
    iso_matrix2 = []
    for s in samp1.split(','):
        name = os.path.basename(s)
        iso_matrix = '{}/ISO_module/{}.IsoMatrix'.format(para.out_dir, name)
        iso_matrix1.append(iso_matrix)
    for s in samp2.split(','):
        name = os.path.basename(s)
        iso_matrix = '{}/ISO_module/{}.IsoMatrix'.format(para.out_dir, name)
        iso_matrix2.append(iso_matrix)
    para.module_iso_matrix1 = ','.join(iso_matrix1)
    para.module_iso_matrix2 = ','.join(iso_matrix2)
    para.module_iso_exon = re.sub(r'Matrix$', 'Exon', iso_matrix1[0])
    return


def EM(para, cur_dir):
    utils.exec_cmd(sys.stderr, 'rMATSEM', 'mkdir {}/EM_out 2> /dev/null'.format(para.out_dir))
    para.EM_out = '{}/EM_out/EM.out'.format(para.out_dir)
    cmd = 'Rscript {}/rmats-EM/IsoRscript.r {} {} {} {} FALSE {} {} {} {} {} {} {}'.format(cur_dir,
                                                                                           para.module_iso_exon,
                                                                                           para.module_iso_matrix1,
                                                                                           para.module_iso_matrix2,
                                                                                           para.EM_out,
                                                                                           para.nthread,
                                                                                           para.em_seed,
                                                                                           para.em_draw_samp,
                                                                                           para.em_p_thres,
                                                                                           para.em_pair_simu,
                                                                                           para.em_n_iter,
                                                                                           para.em_epsilon)
    utils.exec_cmd(sys.stderr, 'rMATSEM', cmd)
    return


def ISO_classifyify(para, cur_dir):
    utils.exec_cmd(sys.stderr, 'ISOClassify', 'mkdir {}/ISO_classify 2> /dev/null'.format(para.out_dir))
    para.module_type_summary = '{}/ISO_classify/ISO_module_type_summary.txt'.format(para.out_dir)
    para.module_type = '{}/ISO_classify/ISO_module_type.txt'.format(para.out_dir)
    cmd = 'python {}/ISOClassify/IsoClass.py {} {} {}'.format(cur_dir, para.module_iso_exon, para.module_type_summary,
                                                              para.module_type)
    utils.exec_cmd(sys.stderr, 'ISOClassify', cmd)
    return


def ISO_plot(para, cur_dir):
    utils.exec_cmd(sys.stderr, 'ISOPlot', 'mkdir {}/ISO_plot 2> /dev/null'.format(para.out_dir))
    para.module_coor = '{}/ISO_classify/ISO_module_coor.txt'.format(para.out_dir)
    para.module_gene = '{}/ISO_classify/ISO_module_gene.txt'.format(para.out_dir)
    cmd = 'python {}/ISOPlot/IsoPlot.py {} {} {} {}'.format(cur_dir, para.EM_out, para.module_iso_exon,
                                                            para.module_coor,
                                                            para.module_gene)
    utils.exec_cmd(sys.stderr, 'IsoPlot', cmd)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="rMATS-ISO: a computational and statistical framework for quantifying mRNA isoform variation using RNA-seq read counts and splicing module information.")
    # "\nquick start: python rMATS-ISO.py --in-gtf test_data/in.gtf --in-bam short_read_bam_input.list.example -o test_out") # TODO example command in README.md
    parse_argv(parser)
    args = parser.parse_args()
    para = iso_para(args)
    cur_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

    # running steps
    run_lr2rmats = True if para.in_long_list and para.in_short_list else False
    if para.in_long_list and (not para.in_short_list or not para.in_genome):
        utils.fatal_format_time('lr2ramts', 'Please provide genome sequence and short-read data to run \'lr2rmats\'.\n')

    run_short_map = False if not para.in_long_list and para.in_bam_list else True

    run_ISO_module = True
    run_em = True
    run_ISO_classifyify = True
    run_iso_plot = True

    utils.exec_cmd(sys.stderr, 'rMATS-ISO', 'mkdir {} 2> /dev/null'.format(para.out_dir))
    # 0. long-read to updated GTF
    # lr2rmats: long-read => GTF
    if run_lr2rmats:
        lr2rmats(para, cur_dir)
        para.in_gtf = para.updated_gtf

    # 1. short-read to BAM
    # STAR: short-read => BAM
    if run_short_map:
        short_map(para)
    # 2. splice module detection and quantification
    # IsoModule: GTF + BAM => .IsoExon/.IsoMatrix
    if run_ISO_module:
        ISO_module(para, cur_dir)
    # 2. statistical test
    # rMATS-EM: .IsoExon/.IsoMatrix => EM.out
    if run_em:
        EM(para, cur_dir)

    # 3. IsoClassify && IsoPlot
    if run_ISO_classifyify:
        ISO_classifyify(para, cur_dir)
    if run_iso_plot:
        ISO_plot(para, cur_dir)
