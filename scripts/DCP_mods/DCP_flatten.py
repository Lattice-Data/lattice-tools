import json
import os
import pandas as pd


def flatten_obj(obj):
	new_obj = {}
	for k,v in obj.items():
		if isinstance(v, dict):
			for x,y in v.items():
				new_obj[k + '.' + x] = y
		elif (isinstance(v, list) and all(isinstance(i, dict) for i in v)):
			new_v = []
			for i in v:
				new_i = flatten_obj(i)
				new_v.append(new_i)
			new_obj[k] = new_v
		elif v in [True, False]:
			new_obj[k] = str(v)
		else:
			new_obj[k] = v
	return new_obj


def get_inputs(links, obj_id):
	obj_id = obj_id.split('/')[-1]
	ins = []
	for l in links['links']:
		if l.get('outputs'):
			for outs in l['outputs']:
				if obj_id in outs['output_id']:
					for inp in l['inputs']:
						ins.append(inp['input_type'] + '/' + inp['input_id'])
					for prot in l['protocols']:
						ins.append(prot['protocol_type'] + '/' + prot['protocol_id'])
	return ins


def tsv_report(ds_id):
	uber_dict = {}
	metadata_dict = {}

	# complile all metadata file contents in a dictionary
	for obj_type in os.listdir(ds_id + '/metadata'):
		for o in os.listdir(ds_id + '/metadata/' + obj_type):
			obj_id = obj_type + '/' + o.split('_')[0]
			metadata_dict[obj_id] = {}
			o_json = json.load(open(ds_id + '/metadata/' + obj_type + '/' + o))
			for k,v in flatten_obj(o_json).items():
				if v:
					if isinstance(v, list) and isinstance(v[0], dict):
						for e in v:
							for k2,v2 in e.items():
								if k + '.' + k2 in metadata_dict[obj_id]:
									if isinstance(v2, list):
										metadata_dict[obj_id][k + '.' + k2].extend(v2)
									else:
										metadata_dict[obj_id][k + '.' + k2].append(str(v2))
								elif isinstance(v2, list):
									metadata_dict[obj_id][k + '.' + k2] = v2
								else:
									metadata_dict[obj_id][k + '.' + k2] = [str(v2)]
					elif isinstance(v, dict):
						for k2,v2 in v.items():
							metadata_dict[obj_id][k + '.' + k2] = v2
					else:
						metadata_dict[obj_id][k] = str(v)

	# cycle through to pull files & walk links backward
	for f in os.listdir(ds_id + '/links'):
		link_json = json.load(open(ds_id + '/links/' + f))

		supps = []
		for l in link_json['links']:
			if l.get('link_type') == 'supplementary_file_link':
				supps = [f['file_type'] + '/' + f['file_id'] for f in l['files']]

		seq_file_link = link_json['links'][0]
		for seq_file in seq_file_link['outputs']:
			if seq_file['output_type'] != 'sequence_file':
				sys.exit('{} in links file {} not sequence_file'.format(seq_file['output_id'], f))
			seq_file_id = seq_file['output_id']
			uber_dict[seq_file_id] = {}
			ins = []
			ins.extend(get_inputs(link_json, seq_file_id))

			seen = set()
			remaining = ins
			while remaining:
				seen.update(remaining)
				next_remaining = set()
				for identifier in remaining:
					ins.extend(get_inputs(link_json, identifier))
				remaining = next_remaining - seen

			for k,v in metadata_dict['sequence_file/' + seq_file_id].items():
				prop = 'sequence_file.' + k
				if prop in uber_dict[seq_file_id]:
					uber_dict[seq_file_id][prop].append(v)
				elif isinstance(v, list):
					uber_dict[seq_file_id][prop] = v
				else:
					uber_dict[seq_file_id][prop] = [v]

			for k,v in metadata_dict['project/' + ds_id].items():
				prop = 'project.' + k
				if prop in uber_dict[seq_file_id]:
					uber_dict[seq_file_id][prop].append(v)
				elif isinstance(v, list):
					uber_dict[seq_file_id][prop] = v
				else:
					uber_dict[seq_file_id][prop] = [v]

			if supps:
				for s in supps:
					for k,v in metadata_dict[s].items():
						prop = 'supplementary_file.' + k
						if prop in uber_dict[seq_file_id]:
							uber_dict[seq_file_id][prop].append(v)
						elif isinstance(v, list):
							uber_dict[seq_file_id][prop] = v
						else:
							uber_dict[seq_file_id][prop] = [v]

			# complile metadata for all input objects up the graph
			for i in ins:
				for k,v in metadata_dict[i].items():
					prop = i.split('/')[0] + '.' + k
					if prop in uber_dict[seq_file_id]:
						if isinstance(v, list):
							v = set(v)
							uber_dict[seq_file_id][prop].append('||'.join(v))
						else:
							uber_dict[seq_file_id][prop].append(v)
					elif isinstance(v, list):
						uber_dict[seq_file_id][prop] = v
					else:
						uber_dict[seq_file_id][prop] = [v]

			for k,v in uber_dict[seq_file_id].items():
				v = set(v)
				uber_dict[seq_file_id][k] = '||'.join(v)

	df = pd.DataFrame(uber_dict).fillna('').transpose()
	df.to_csv('DCP_outs/' + ds_id + '.tsv', sep='\t')
