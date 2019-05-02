import re

import pandas as pd
from collections import defaultdict
standard_data = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/RS_judge/resistant_standard.tab", sep='\t', header=None)
standard = standard_data.iloc[3:, 6:-2]
standard.columns = ['name1',
                    'name2',
                    'disk content',
                    'S',
                    'I',
                    'R',
                    'S.1',
                    'I.1',
                    'R.1',
                    ]


def generate_fun(val):
    if '≥' in val:
        val = float(val.split('≥')[-1])
        return (lambda x: x >= val, val)
    elif '≤' in val:
        val = float(val.split("≤")[-1])
        return (lambda x: x <= val, val)
    elif '-' in val:
        _min, _max = sorted(map(float, val.split('-')))
        if _min == _max:
            return (lambda x: x == _min, _min)
        else:
            return (lambda x: x >= _min and x <= _max, (_min,_max))
    elif re.match('^[0-9]+$', val):
        return (lambda x: x == float(val), val)
    else:
        import pdb;pdb.set_trace()
        raise Exception


def process_SR(col):
    standard = {}
    for t in col.index:
        val = col[t]
        if type(val) == str and str(val).strip() != '-':
            val = str(val).strip()
            if '/' in val:
                get_num = re.findall("[0-9]+", val)
                if '-' in val:
                    finalvals = ["%s-%s" % (get_num[0], get_num[2]),
                                 "%s-%s" % (get_num[1], get_num[3])]

                else:
                    finalvals = [val[0] + _ for _ in get_num]
                standard[t] = [generate_fun(_) for _ in finalvals]
            else:
                standard[t] = generate_fun(val)
    return standard


standard_dict = defaultdict(dict)
for idx, row in standard.iterrows():
    for name in row.loc[['name1', 'name2']]:
        if type(name) == str and bool(str(name)):
            name = name.strip()
            format_names = re.findall('[a-zA-Z]+', name)
            for n in format_names:
                standard_dict[n]['disk'] = process_SR(row.loc[['S',
                                                               'I',
                                                               'R']])
                standard_dict[n]['MIC'] = process_SR(row.loc[['S.1',
                                                              'I.1',
                                                              'R.1']])

    ############################################################

ori_data = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/RS_judge/total_actineto.csv", sep='\t', index_col=1)

ori_data = ori_data.loc[~ori_data.index.isnull(), :]
ori_data = ori_data.loc[~ori_data.index.duplicated(keep=False), :]

med = ['AMC', 'AMK', 'AMK.1', 'AMP', 'ATM', 'CAZ', 'CAZ.1', 'CEC', 'CIP',
       'CIP.1', 'CRO', 'CRO.1', 'CSL', 'CSL.1', 'CTX', 'CZO', 'FEP',
       'FEP.1', 'FOX', 'GEN', 'GEN.1', 'IPM', 'IPM.1', 'LVX', 'MEM',
       'MEM.1', 'MNO', 'MNO.1', 'NAL', 'NIT', 'PIP', 'POL', 'POL.1',
       'SAM', 'SAM.1', 'SXT', 'SXT.1', 'TGC', 'TGC.1', 'TZP', 'TZP.1']


############################################################
remained_med = []
for i in med:
    if i.replace('.1','') in standard_dict:
        remained_med.append(i)
med_data = ori_data.loc[:,remained_med]


for idx,val in med_data.iteritems():
    if '.1' not in idx:
        j_funcs = standard_dict[idx]['disk']
    else:
        j_funcs = standard_dict[idx.replace('.1','')]["MIC"]
    status = ['no'] * len(val)
    new_val = []
    for i in val:
        if type(i) ==str:
            new_val.append(float(re.findall("(0?\.?[0-9]+)",i)[0]))
        else:
            new_val.append(float(i))
    for nid,v in enumerate(new_val):

        for rs in j_funcs.keys():
            j_func = j_funcs[rs][0]
            if type(j_func) == list:
                continue
            if type(j_func) == tuple:
                j_func = j_func[0]
            if j_func(v):
                status[nid] = rs.replace('.1','')
    med_data.loc[:,"RES_%s" % idx] = status

############################################################
# filter and merge
dropped_iterms = []
merged_col = defaultdict(list)
for idx,val in med_data.iteritems():
    if idx.startswith('RES_'):
        if (val == 'no').all():
            dropped_iterms.append(idx)
        else:
            merged_col[idx.replace('.1','')].append(idx)
for k,v in merged_col.items():
    if len(v) == 2:
        c1 = med_data.loc[:,v[0]]
        c2 = med_data.loc[:,v[1]]
        new_v = []
        for v1,v2 in zip(c1,c2):
            if v1 ==v2:
                new_v.append(v1)
            elif v1 != 'no' and v2 != 'no' and v1 !=v2:
                new_v.append('%s;%s' % (v1, v2))
            elif v1 == 'no' or v2 == 'no':
                _v = {'R','S','I','no'}.intersection({v1,v2}).difference({'no'})
                if _v:
                    new_v.append(_v.pop())
                else:
                    new_v.append('no')
            else:
                raise Exception
        dropped_iterms.append(set(v).difference({k}).pop())
        med_data.loc[:,k] = new_v

med_data = med_data.drop(dropped_iterms,axis=1)
med_data.to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/RS_judge/RS.csv",sep=',')






