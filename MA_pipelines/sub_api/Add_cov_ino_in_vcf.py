import pandas as pd
import vcf
from tqdm import tqdm


def parsed_vcf2pos_list(vcf_path):
    """
    Read a vcf file and turn it into a 0-coordinated pos list.
    :param vcf: vcf file path
    :return: a 0-coordinated pos list.
    """
    if type(vcf_path) == str:
        vcf_readed = vcf.Reader(open(vcf_path, 'r'))
    else:
        try:
            vcf_readed = vcf.Reader(fsock=vcf_path)
        except:
            raise IOError('Wrong vcf, it is a %s' % str(type(vcf)))

    pos_list = []
    for record in vcf_readed:
        if record.is_snp:
            pos_list.append((record.CHROM,
                             record.POS - 1))
    return pos_list


def count_sample(samples):
    if '-' in samples[0]:
        # todo: may got wrong. using '-' as separator to distinguish samples
        return list(set([sample.rpartition('-')[0]
                         for sample in samples]))
    else:
        return list(set(samples))


def Add_in_vcf_SO(infofile,
                  vcf_path,
                  output_vcf, ):
    """
    receive a vcf file and a related bam.
    Add coverage from bam into vcf and make it a new field.
    For Tumor/Normal only or germline vcf.

    :param bam: bam path
    :param vcf: vcf path or vcf file object.
    :return: A new vcf to output.
    """
    ori_format2info = ['AF', 'AD']

    if type(vcf_path) == str:
        vcf_readed = vcf.Reader(open(vcf_path, 'r'))
    else:
        try:
            vcf_readed = vcf.Reader(fsock=vcf_path)
        except:
            raise IOError('Wrong vcf, it is a %s' % str(type(vcf_path)))

    info_df = pd.read_csv(infofile, sep='\t')
    info_df.index = info_df.iloc[:, [1, 2, 3]].astype(str).sum(1)
    samples = count_sample(vcf_readed.samples)

    if len(samples) != 1:
        return

    new_infos = vcf_readed.infos
    machine = list(new_infos.values())[0]

    field1 = "SAD"
    field2 = "SAF"
    field3 = "PoS"
    field1_info = [field1, 'R', 'Integer',
                   "(REF base count, alt base count). Self cal coverage from bam file. It is different from AD. It is rawer than AD."]
    field2_info = [field2, 'R', 'Float',
                   "Alt base count divide the total reads number in this pos. Self cal frequency from bam file. It is different from AF. It is rawer than AF."]
    field3_info = [field3, '.', 'Integer',
                   "A field which describe this file is single only analysis or pair analysis. 1 for single analysis, 2 for pair analysis."]
    new_infos[field1] = machine._make(field1_info + [None, None])
    new_infos[field2] = machine._make(field2_info + [None, None])
    new_infos[field3] = machine._make(field3_info + [None, None])
    for ori_format in ori_format2info:
        if ori_format not in new_infos.keys():
            new_infos[ori_format] = list(new_infos.values())[0]._make(
                list(vcf_readed.formats[ori_format]._asdict().values()) + [None, None])

    vcf_readed.infos = new_infos

    vcf_writer = vcf.Writer(open(output_vcf, 'w'),
                            vcf_readed)

    for record in tqdm(vcf_readed):
        if record.is_snp:
            # SNP instead of indel
            query_index = record.CHROM + str(record.POS - 1) + str(record.REF)
            if query_index not in info_df.index:
                # it means there are np reads here
                ref_cov, alt_cov = 0, 0
            else:
                row = info_df.loc[query_index, :]
                if len(row.shape) == 2:
                    # if multiple index occur
                    row = row.iloc[0, :]

                ref_base = row["Reference"]  # should same as record.REF
                ref_cov = row[ref_base.upper()]
                alt_cov = row[str(record.ALT[0]).upper()]

            SAD = [int(ref_cov),
                   int(alt_cov)]
            try:
                SAF = round(float(alt_cov) / sum(SAD), 4)
            except ZeroDivisionError:
                SAF = 0

            record.INFO[field1] = tuple(SAD)
            record.INFO[field2] = SAF
            record.INFO[field3] = 1
        else:
            # for indel we just ignore it.
            # write the original info
            pass

        for sample in record.samples:
            data = dict(sample.data._asdict())
            for ori_format in ori_format2info:
                if data.get(ori_format,''):
                    record.INFO[ori_format] = data[ori_format]

        vcf_writer.write_record(record)
    vcf_writer.close()

#
# def Add_in_vcf_PA(bam_list, vcf_path, output_vcf,
#                   fasta_file='/home/liaoth/data/hg19/ucsc.hg19.fasta'):
#     """
#     receive a vcf file and a related bam.
#     Add coverage from bam into vcf and make it a new field.
#     For pair analyse vcf.
#     bam_list order must like [normal one, tumor one]
#     :param bam: bam path
#     :param vcf: vcf path or vcf file object.
#     :return: A new vcf to output.
#     """
#     ori_format2info = ['AF','AD']
#     field1 = "SAD"
#     field2 = "SAF"
#     field3 = "PoS"
#     NT_SIG = [pfn(_bam, 'mt2_for') for _bam in bam_list]
#     NT_name = [pfn(_bam, 'sample_name') for _bam in bam_list]
#     if type(vcf_path) == str:
#         vcf_readed = vcf.Reader(open(vcf_path, 'r'))
#     else:
#         try:
#             vcf_readed = vcf.Reader(fsock=vcf_path)
#         except:
#             raise IOError('Wrong vcf, it is a %s' % str(type(vcf_path)))
#
#     pos_list = parsed_vcf2pos_list(vcf_path)
#
#     is_single = False
#     right_infos = vcf_readed.infos
#     machine = right_infos.values()[0]
#     # Modify the info part.
#     if not is_single:
#         field1_info = [field1, '4', 'Integer',
#                            "(REF base count, alt base count). Self cal allele depths from bam file. If there are two pair, it is normal-tumore order."]
#         field2_info = [field2, 'R', 'Float',
#                            "Alt base count divide the total reads number in this pos. Self cal frequency from bam file. If there are two pair, it is normal-tumore order."]
#         field3_info = [field3, '1', 'Integer',
#                        "A field which describe this file is single only analysis or pair analysis. 1 for single analysis, 2 for pair analysis."]
#
#         right_infos[field1] = machine._make(field1_info + [None, None])
#         right_infos[field2] = machine._make(field2_info + [None, None])
#         right_infos[field3] = machine._make(field3_info + [None, None])
#         for ori_format in ori_format2info:
#             _ori_format_info = vcf_readed.formats[ori_format]._asdict().values() + [None,None]
#             _ori_format_info[3]+= ". If there are two pair, it is normal-tumore order." # fetch ori format value and ID and fix it into length == 6
#             if ori_format == 'AD':
#                 _ori_format_info[1] = '4'
#             elif ori_format == 'AF':
#                 _ori_format_info[1] = 'R'
#             right_infos[ori_format] = machine._make(_ori_format_info)
#
#     vcf_readed.infos = right_infos
#     # Fetch the cov info from bam file, and prepare the writed file.
#     all_cov_info = special_cal_cov(bam_list, pos_list, fasta_file)
#     vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_readed)
#
#     for record in vcf_readed:
#         if record.is_snp:
#             query_for = (record.CHROM, record.POS - 1)
#             buckec_SAD = []
#             bucket_SAF = []
#             for ori_format in ori_format2info:
#                 exec('bucket_%s = []' % ori_format)
#
#             for sample_call in record.samples:
#                 # it needs to fix the sample and the cov_info order.
#                 sample = str(sample_call.sample)
#
#                 idx = [NT_SIG.index(s) for s, n in zip(NT_SIG, NT_name) if sample == n][0]
#                 cov_info = all_cov_info[idx][query_for]
#
#                 ref_base, ref_cov = cov_info[0]
#                 if len(cov_info) > 2:
#                     for n_i in range(1, len(cov_info) - 1):
#                         if cov_info[n_i][0] == record.ALT[0]:
#                             alt_base, alt_cov = cov_info[n_i]
#                 elif len(cov_info) == 1:
#                     alt_base = record.ALT[0]
#                     alt_cov = 0
#                 else:
#                     alt_base, alt_cov = cov_info[1]
#                 ### fix the bucket order to normal-tumore order.
#                 if sample == [n for s, n in zip(NT_SIG, NT_name) if s == N_sig]:
#                     buckec_SAD.insert(0,int(alt_cov))
#                     buckec_SAD.insert(0,int(ref_cov))
#
#                     if sum((int(ref_cov), int(alt_cov))) != 0:
#                         bucket_SAF.insert(0,round(float(alt_cov) / sum((int(ref_cov), int(alt_cov))), 4))
#                     else:
#                         bucket_SAF.insert(0,0)
#                     # data = dict(sample_call.data._asdict())
#                     for ori_format in ori_format2info:
#                         if ori_format =='AD':
#                             exec("bucket_{i}.insert(0,tuple(data['{i}'])[0])".format(i=ori_format))
#                             exec("bucket_{i}.insert(0,tuple(data['{i}'])[1])".format(i=ori_format))
#                         else:
#                             exec("bucket_{i}.insert(0,data['{i}'])".format(i=ori_format))
#                 else:
#                     buckec_SAD+=[int(ref_cov), int(alt_cov)]
#
#                     if sum((int(ref_cov), int(alt_cov))) != 0:
#                         bucket_SAF.append(round(float(alt_cov) / sum((int(ref_cov), int(alt_cov))), 4))
#                     else:
#                         bucket_SAF.append(0)
#                     # data = dict(sample_call.data._asdict())
#                     for ori_format in ori_format2info:
#                         if ori_format == 'AD':
#                             exec("bucket_{i} += list(data['{i}'])".format(i=ori_format))
#                         else:
#                             exec("bucket_{i}.append(data['{i}'])".format(i=ori_format))
#             record.INFO[field1] = buckec_SAD
#             record.INFO[field2] = bucket_SAF
#             record.INFO[field3] = 2
#             for ori_format in ori_format2info:
#                 exec("record.INFO['{i}'] = bucket_{i}".format(i=ori_format))
#         vcf_writer.write_record(record)
#     vcf_writer.close()
