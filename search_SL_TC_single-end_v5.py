import sys
import re
import subprocess
from collections import defaultdict
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import time

# Setup the program arguments
program = 'search_SL_TC_single-end_v5.py'
parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description='''Search for SL sequences in single-end fastq reads and filter potential hits...

	V5 has the final program name, sample_id required for temporal files and reduce redundancy when running multiple jobs on the same path
    Run multithread mode in Linux, MAC and Windows (Tested in Mac-intel and Linux)
    Dependencies: Bowtie2, gffread (need to be installed and available in the path)
    
    Tested on python 3.8.12
    ''')

requiredNamed = parser.add_argument_group('Mandatory arguments')
requiredNamed.add_argument("-s", "--spliced", dest='sl_seq', required=True, type=str, help='The reference SL sequence in fasta format')
requiredNamed.add_argument("-q", "--fastq", dest='fqfile', required=True, type=str, help='Fastq reads, single end')
requiredNamed.add_argument("-g", "--genome", dest='refgenom', required=True, type=str, help='Genome reference, fasta format')
requiredNamed.add_argument("-i", "--gindex", dest='g_index', required=True, type=str, help='Prefix for the Bowtie genome index')
optionalNamed = parser.add_argument_group('Optional arguments')
optionalNamed.add_argument("-t", "--threads", dest='threads', default=1, type=int, help='Number of threads to use, default=1')
optionalNamed.add_argument("-p", "--sprefix", dest='sample_id', default='search_sl_trimmer', type=str, help='Sample ID prefix, default=search_sl_trimmer. Use this parameter if youre running multiple jobs on the same path, the script creates multiple temporary files that will cause redundant files')
optionalNamed.add_argument("-m", "--mismatch", dest='mmatchs', default=6, type=int, help='Maximum number of mismatches in the SL, default=6')
optionalNamed.add_argument("-l", "--slength", dest='slength', default=8, type=int, help='Minimum length of the SL, default=8')

args = parser.parse_args()

# Load SL sequence
with open(args.sl_seq, 'r') as leader_file:
    leader_name = leader_file.readline()[1:].rstrip('\n')
    leader = leader_file.readline().rstrip('\n').upper()

def refasta(infile):
	""" Reformatea un fichero fasta para que la secuencia esté en una sóla línea """
	with open(infile, 'r') as file:
		s = ''
		seq = {}
		for line in file:
			line = line.rstrip('\n')
			if line:
				if line.startswith('>'):
					if s:
						seq[name] = s
						s = ''
					name = line
				else:
					s = s + line
		seq[name] = s

	with open (infile, 'w') as file:
		for s in seq.keys():
			file.write(s + '\n')
			file.write(seq[s] + '\n')

# Define constants
min_leader_length = args.slength
max_mm = args.mmatchs
max_error = max_mm / len(leader)
bow = 'bowtie2 -p ' + str(args.threads) + ' --np 0 --n-ceil L,0,0.02 --rdg 0,6 --rfg 0,6 --mp 6,2 --score-min L,0,-0.24'
genome = args.refgenom

# Thread-safe dictionary to store potential SL hits
sl = defaultdict(list)

def align(query, subject):
    """ Compares two sequences and calculates the mismatch percentage """
    mm = sum(1 for q, s in zip(query, subject) if q != s)
    return mm / len(query)

def process_fastq_chunk(chunk, leader, max_error):
    """ Process a chunk of the fastq file to find SL hits """
    potential_hits = []
    local_sl = defaultdict(list)

    for record in chunk:
        name, seq, coment, qual = record

        for i in range(0, len(leader) - min_leader_length + 1):
            query = leader[i:]
            subject = seq[:len(leader) - i]

            # If align score is within max error, consider it a potential hit
            if align(list(query), list(subject)) <= max_error:
                local_sl[name[1:].split(' ')[0]] = [name, seq, coment, qual, len(query), 0]
                potential_hits.append((name, seq[len(query):], coment, qual[len(query):]))
                break

    return potential_hits, local_sl

def read_fastq_in_chunks(fqfile, chunk_size=100000):
    with open(fqfile, 'r') as fastq_infile:
        chunk = []
        while True:
            name = fastq_infile.readline().rstrip('\n')
            seq = fastq_infile.readline().rstrip('\n')
            coment = fastq_infile.readline().rstrip('\n')
            qual = fastq_infile.readline().rstrip('\n')

            if not name:
                if chunk:
                    yield chunk
                break

            chunk.append((name, seq, coment, qual))
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []

def main():
	print('\n\nSL search in fastq file has started...')
	start_time = time.perf_counter()
    
	with ProcessPoolExecutor(max_workers=args.threads) as executor:
		futures = []
		for chunk in read_fastq_in_chunks(args.fqfile):
			futures.append(executor.submit(process_fastq_chunk, chunk, leader, max_error))

		with open(args.sample_id + '_seqleader_temp.fastq', 'w') as leader_outfile:
			count = 0
			for future in as_completed(futures):
				potential_hits, local_sl = future.result()
				count += len(potential_hits)

				for record in potential_hits:
					leader_outfile.write('\n'.join(record) + '\n')

				sl.update(local_sl)
    
	end_time = time.perf_counter()
	elapsed_time = (end_time - start_time) / 60
	print(f"Time elapsed: {elapsed_time:.2f} minutes")
	print(f"Potential Leader hits: {count}" + '\n')

	print('Bowtie2 alignment has started has started...\n')
	subprocess.call(f"{bow} -x {args.g_index} -U " + args.sample_id + "_seqleader_temp.fastq -S " + args.sample_id + "_temp.sam", shell=True)
	print('\nBowtie2 alignment has finished...\n')

	##########################################################################################

	#--Recorremos el fichero sam para extraer las secuencias genómicas inmediatamente anteriores al inicio del alineamiento de los potenciales SL
	#--De este modo podemos analizar si la secuencia leader estaba en el genoma y por lo tanto es un falso positivo
	print('SAM analisis has started...')
	with open(args.sample_id + '_temp.gtf', 'w') as gtf_temp, open(args.sample_id + '_temp.sam', 'r') as sam:
		chr = {}

		for line in sam:
			line = line.rstrip('\n')
			col = line.rstrip('\n').split('\t')
		
			#--skip header lines and take chromosome sizes
	
			if '@' in line[0]: 
				if col[0] == '@SQ':
					chr[col[1].split(':')[1]] = int(col[2].split(':')[1])
					continue
				else:
					continue

			#--The read do not align
			if col[2] == '*':		
				del sl[col[0]]
				continue
		
			#--The read align
			else: 
				#--Obtenermos la información del bitwise flag
				bitwise = int(col[1])
				bits = bin(bitwise)[2:]
				byte = "{0:0>12}".format(bits)

				#--Obtenemos las coordenadas en función de la orientación del alineamiento
				if int(byte[-5]): #--int(byte[-5])==1 (reverse)
					strand = '-'					
					cigar = col[5]
					let = re.split('\d+', cigar)[1:] #--Dividimos por número	
					val = re.split('\D+', cigar)[:-1] #--Dividimos por letras

					coord = 0
					for l, v in zip(let, val):
						if l == 'M': #--Fragmento de lectura alineado
							coord = coord + int(v)
						elif l == 'D': #--Delección en la lectura
							coord = coord + int(v)
						else: #--Las inserciones no las sumamos
							continue
					coord = int(col[3]) + coord - 1
			
					x = coord + 1
					y = coord + sl[col[0]][-2]
			
					if y > chr[col[2]]:
						y = chr[col[2]]
						
					if (y - x) + 1 < min_leader_length: #--Si la secuencia esta fuera de los límites de cromosoma damos por bueno al SL, porque no podemos discriminarlo
						sl[col[0]][-1] = 1
			
					else:		
						#--Escribimos las coordenadas del fragmento eliminado en un fichero gtf para extraer la secuencia genómica
						transcript = [col[2], 'CBMSO', 'transcript', str(x), str(y), '.', '-', '.']
						att = ['gene_id "', col[0], '"; ', 'transcript_id "',col[0], '"; '] 	
						gtf_temp.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
	
				else: #--int(byte[-5])==0 (forward)
					strand = '+'
					coord = int(col[3])				
			
					y = coord - 1
					x = coord - sl[col[0]][-2]
			
					if x <= 0:
						x = 1
			
					if (y - x) + 1 < min_leader_length:  #--Si la secuencia esta fuera de los límites de cromosoma damos por bueno al SL, porque no podemos discriminarlo
						sl[col[0]][-1] = 1
				
					else:
						#--Escribimos las coordenadas del fragmento eliminado en un fichero gtf para extraer la secuencia genómica
						transcript = [col[2], 'CBMSO', 'transcript', str(x), str(y), '.', '+', '.']
						att = ['gene_id "', col[0], '"; ', 'transcript_id "',col[0], '"; '] 	
						gtf_temp.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')

	print('SAM analisis has finished...\n')	
	#--Extraemos las secuencias genómicas
	subprocess.call('gffread -w ' + args.sample_id + '_temp.fa -g ' + genome + ' ' + args.sample_id + '_temp.gtf', shell = True)

	#--Reformateo del fichero temp.fa por si la secuencia ocupa de más de 1 línea (en el caso de que la secuencia eliminada sea mayor de 60 pb)
	refasta(args.sample_id + '_temp.fa')
	print('Search of SL false positives has started...')
	#--Analizamos las secuencias genómicas en busca de secuencias leader genómicas (mismas condiciones que en las lecturas)
	with open (args.sample_id + '_temp.fa', 'r') as seqs:
		while True:
			name = seqs.readline().rstrip('\n')[1:]
			seq = seqs.readline().rstrip('\n')
 	
			if not name: #--Finish
				break
 	
			query = leader[-(sl[name][-2]):]
		
			if len(seq) < len(query): #--para las secuencias de los extremos en los que falta secuencia genómica, aunque siempre hay más que el mínimo (filtrado arriba)
				query = leader[-(len(seq)):]
				if align(list(query), list(seq)) > max_error:
					sl[name][-1] = 1 #--Cambiamos la etiqueta del diccionaria a 1, es decir, que es un SL genuino y debe ser recortado	
			else:
				if align(list(query), list(seq)) > max_error:
					sl[name][-1] = 1 #--Cambiamos la etiqueta del diccionaria a 1, es decir, que es un SL genuino y debe ser recortado

	#--Volvemos a recorrer el fichero con las lecturas originales y ahora si recortamos los SL que sabemos son correctos
	sl_seqs_file = args.sample_id + '_SL_seqs.fasta'
	fastq_trimmed = args.sample_id + '_SL_trimmed.fastq'
	
	with open (args.fqfile, 'r') as fastq_infile, open(fastq_trimmed, 'w') as reads_outfile, open(sl_seqs_file, 'w') as out_sl: #open('seqleader_ids', 'w') as leader_outfile, 
		count = 0
		while True:
			name = fastq_infile.readline().rstrip('\n')
			seq = fastq_infile.readline().rstrip('\n')
			coment = fastq_infile.readline().rstrip('\n')
			qual = fastq_infile.readline().rstrip('\n')
	
			if not name: #--Finish
				break

			if sl[name[1:].split(' ')[0]] and sl[name[1:].split(' ')[0]][-1]:
				count = count + 1
				new_name = name.split(' ')[0] + 'L ' + ' '.join(name.split(' ')[1:])
				#new_name = name + 'L' #--Añado una L para identificar las lecturas después
				record = [new_name, seq[sl[name[1:].split(' ')[0]][-2]:], coment, qual[sl[name[1:].split(' ')[0]][-2]:]] 
				#leader_outfile.write(name[1:].split(' ')[0] + '\n') #--Guardar sólo los ids y luego sacarlos del sam general
				reads_outfile.write('\n'.join(record) + '\n')
				out_sl.write('>' + name[1:] + '\n' + seq[:sl[name[1:].split(' ')[0]][-2]] + '\n')
			else:
				record = [name, seq, coment, qual] 
				reads_outfile.write('\n'.join(record) + '\n')

	print('Search of SL false positives has finished...\n')
	print(f'File containing SL sequences for MSA has been created: {sl_seqs_file}')
	print(f'File containing fastq reads with SL trimmed has been created: {fastq_trimmed}')	
	print(f'Valid leader hits: {count}')
	print('\n\n')
	os.remove(args.sample_id + '_seqleader_temp.fastq')
	os.remove(args.sample_id + '_temp.sam')
	os.remove(args.sample_id + '_temp.gtf')
	os.remove(args.sample_id + '_temp.fa')

if __name__ == '__main__':
    main()
