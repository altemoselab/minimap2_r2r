import sys
import numpy as np 

#m64076_211222_124721/<integer>/1_10000 for subread
#m64076_211222_124721/<integer>/ccs for ccs



ref = sys.argv[1]
query = sys.argv[2]

conversion_dict = {}

ref_out = ref[:-2] + "map_compatible.fa"
ref_out_file = open(ref_out,'w')

# randomly select 17 digits for smrt cell m+5+6+6
digits = np.random.choice(range(0,10),size=17,replace=True).astype(str)
smrt_cell = "m" + ''.join(digits[:5]) + "_" + ''.join(digits[5:11]) + "_" + ''.join(digits[11:])


count = 1
with open(ref,'r') as handle:
	for line in handle:
		if ">" in line:
			conversion_dict[line.split()[0].strip()[1:]] = count
			ccs_header = '{}/{}/ccs'.format(smrt_cell,str(count))
			ref_out_file.write(">"+ccs_header+"\n")
			count += 1
		else:
			ref_out_file.write(line)

ref_out_file.close()


query_out = query[:-5] + "map_compatible.fastq"
query_out_file = open(query_out,'w')


subread_counts = {}

with open(query,'r') as handle:
	for line in handle:
		if line[0] == "@" and len(line) < 200:
			zmw_label = conversion_dict[line.split()[0].strip()[1:]]
			if zmw_label in subread_counts:
				subread_counts[zmw_label] +=1
			else:
				subread_counts[zmw_label] = 1

			sub_spec_start = subread_counts[zmw_label]
			sub_spec_end = sub_spec_start + 1
			sub_spec_start, sub_spec_end = str(sub_spec_start), str(sub_spec_end)

			subread_header = '{}/{}/{}_{}'.format(smrt_cell,str(zmw_label),sub_spec_start,sub_spec_end)


			query_out_file.write("@"+subread_header+"\n")
			count += 1
		else:
			query_out_file.write(line)


query_out_file.close()


for key in conversion_dict:
	print(key,"\t",smrt_cell,conversion_dict[key])



