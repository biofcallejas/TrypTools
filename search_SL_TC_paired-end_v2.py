import sys
import re
import subprocess
from collections import defaultdict
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import time

# Argument setup
program = 'search_SL_TC_paired-end_v2.py'
parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter,
                                 description='''Search for SL sequences in paired-end fastq reads and filter potential hits...
	Dependencies: Bowtie2, gffread (need to be installed and available in the path)
	
	v2 runs in multithread mode
	Tested on Python 3.8.12
''')

# Define arguments
requiredNamed = parser.add_argument_group('Mandatory arguments')
requiredNamed.add_argument("-s", "--spliced", dest='sl_seq', required=True, type=str, help='The reference SL sequence in fasta format')
requiredNamed.add_argument("-r1", "--fastq1", dest='fqfile_r1', required=True, type=str, help='Fastq reads mate1; fragment r1')
requiredNamed.add_argument("-r2", "--fastq2", dest='fqfile_r2', required=True, type=str, help='Fastq reads mate2; fragment r2')
requiredNamed.add_argument("-g", "--genome", dest='refgenom', required=True, type=str, help='Genome reference, fasta format')
requiredNamed.add_argument("-i", "--gindex", dest='g_index', required=True, type=str, help='Prefix for the Bowtie genome index')
optionalNamed = parser.add_argument_group('Optional arguments')
optionalNamed.add_argument("-t", "--threads", dest='threads', default=1, type=int, help='Number of threads to use, default=1')
optionalNamed.add_argument("-m", "--mismatch", dest='mmatchs', default=6, type=int, help='Maximum number of mismatches in the SL, default=6')
optionalNamed.add_argument("-p", "--sprefix", dest='sample_id', default='search_sl_trimmer', type=str, help='Sample ID prefix, default=search_sl_trimmer. Use this parameter if you are running multiple jobs on the same path, the script creates multiple temporary files that will cause redundant files')
optionalNamed.add_argument("-l", "--slength", dest='slength', default=8, type=int, help='Minimum length of the SL, default=8')
args = parser.parse_args()

def align(query, subject):
    mm = sum(1 for q, s in zip(query, subject) if q != s)
    return mm / len(query)
# Load SL sequence
with open(args.sl_seq, 'r') as leader_file:
    leader_name = leader_file.readline()[1:].rstrip('\n')
    leader = leader_file.readline().rstrip('\n').upper()

min_leader_length = args.slength
max_mm = args.mmatchs
max_error = max_mm / len(leader)
bow = 'bowtie2 -p ' + str(args.threads) + ' --np 0 --n-ceil L,0,0.02 --rdg 0,6 --rfg 0,6 --mp 6,2 --score-min L,0,-0.24'
genome = args.refgenom


def refasta(infile):
    with open(infile, 'r') as file:
        s = ''
        seq = {}
        for line in file:
            line = line.rstrip('\n')
            if line:
                if line[0] == '>':
                    if s:
                        seq[name] = s
                        s = ''
                    name = line
                else:
                    s += line
        seq[name] = s
    with open(infile, 'w') as file:
        for s in seq.keys():
            file.write(s + '\n')
            file.write(seq[s] + '\n')

# Thread-safe dictionary to store potential SL hits
sl = defaultdict(list)

def process_chunk(chunk, leader, min_leader_length, max_error):
    local_sl = defaultdict(list)
    count = 0
    processed_records = []

    for name1, seq1, coment1, qual1, name2, seq2, coment2, qual2 in chunk:
        rl = len(seq1)
        trim = 0
        for i in range(0, len(leader) - min_leader_length + 1):
            query = leader[i:]
            subject = seq1[:len(leader) - i]

            if align(list(query), list(subject)) <= max_error:
                count += 1
                record1 = [name1, seq1[len(query):], coment1, qual1[len(query):]]
                record2 = [name2, seq2, coment2, qual2]
                local_sl[name1.split(' ')[0][1:]] = [name1, seq1, coment1, qual1, 1, len(query), 0]
                processed_records.append((record1, record2))
                trim = 1
                break

        if not trim:
            for i in range(0, len(leader) - min_leader_length + 1):
                query = leader[i:]
                subject = seq2[:len(leader) - i]

                if align(list(query), list(subject)) <= max_error:
                    count += 1
                    record1 = [name1, seq1, coment1, qual1]
                    record2 = [name2, seq2[len(query):], coment2, qual2[len(query):]]
                    local_sl[name2.split(' ')[0][1:]] = [name2, seq2, coment2, qual2, 2, len(query), 0]
                    processed_records.append((record1, record2))
                    trim = 1
                    break
    return processed_records, local_sl

def read_fastq_in_chunks(fqfile_r1, fqfile_r2, chunk_size=100000):
    with open(fqfile_r1, 'r') as fastq_infile1, open(fqfile_r2, 'r') as fastq_infile2:
        chunk = []
        while True:
            name1 = fastq_infile1.readline().rstrip('\n')
            seq1 = fastq_infile1.readline().rstrip('\n')
            coment1 = fastq_infile1.readline().rstrip('\n')
            qual1 = fastq_infile1.readline().rstrip('\n')

            name2 = fastq_infile2.readline().rstrip('\n')
            seq2 = fastq_infile2.readline().rstrip('\n')
            coment2 = fastq_infile2.readline().rstrip('\n')
            qual2 = fastq_infile2.readline().rstrip('\n')

            if not name1:
                break

            chunk.append((name1, seq1, coment1, qual1, name2, seq2, coment2, qual2))

            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []

        if chunk:
            yield chunk

def main():
	print('\n\nSL search in fastq paired-end files has started...')
	start_time = time.perf_counter()
    
	with ProcessPoolExecutor(max_workers=args.threads) as executor:
		futures = []
		for chunk in read_fastq_in_chunks(args.fqfile_r1, args.fqfile_r2):
			futures.append(executor.submit(process_chunk, chunk, leader, min_leader_length, max_error))

	with open(args.fqfile_r1, 'r') as fastq_infile1, open (args.fqfile_r2, 'r') as fastq_infile2, open(args.sample_id + '_seqleader_temp_1.fastq', 'w') as leader_outfile1, open(args.sample_id + '_seqleader_temp_2.fastq', 'w') as leader_outfile2:
		count = 0
		for future in as_completed(futures):
			processed_records, local_sl = future.result()
			count += len(processed_records)
			
			for record in processed_records:
				leader_outfile1.write('\n'.join(record[0]) + '\n')
				leader_outfile2.write('\n'.join(record[1]) + '\n')
			sl.update(local_sl)

	end_time = time.perf_counter()
	elapsed_time = (end_time - start_time) / 60
	print(f"Time elapsed: {elapsed_time:.2f} minutes")
	print(f"Potential Leader hits: {count}" + '\n')
	
##########################################################################################

	# Alignment with Bowtie2 using the dynamic number of threads.
	print('Bowtie2 alignment has started has started...\n')
	#--Alineamos los potenciales SL detectados
	subprocess.call(bow + ' -x ' + args.g_index + ' -1 ' + args.sample_id + '_seqleader_temp_1.fastq -2 ' + args.sample_id + '_seqleader_temp_2.fastq -S ' + args.sample_id + '_temp.sam', shell=True)
	print('\nBowtie2 alignment has finished...\n')

	##########################################################################################

	#--Recorremos el fichero sam para extraer las secuencias genómicas inmediatamente anteriores al inicio del alineamiento de los potenciales SL
	#--De este modo podemos analizar si la secuencia leader estaba en el genoma y por lo tanto es un falso positivo
	print('SAM analisis has started...')
	with open(args.sample_id + '_temp.gtf', 'w') as gtf_temp, open (args.sample_id + '_temp.sam', 'r') as sam:
		chr = {}

		#--leemos la cabecera
		while True:
			line = sam.readline()
			line = line.rstrip('\n')
			col = line.rstrip('\n').split('\t')	
	
			if col[0] == '@SQ':
				chr[col[1].split(':')[1]] = int(col[2].split(':')[1])
				continue
	
			if col[0] == '@PG': #--Para leer sólo la cabecera
				break

		#--Leemos la parte del alineamiento
		while True:	
			line1 = sam.readline().rstrip('\n')
			line2 = sam.readline().rstrip('\n')
	
			if not line1: #--Final del fichero
				break
	
			read_name = line1.split('\t')[0]

			#--Analizamos el bitwise flag de la primera lectura que aparece y determinamos cual de las dos es el pair 1 y 2
			col = line1.split('\t')
			bitwise = int(col[1])
			bits = bin(bitwise)[2:]
			byte = "{0:0>12}".format(bits)
		
			if int(byte[-7]) == 1: #--Par 1
				pair1 = line1
				pair2 = line2	
			else:
				pair1 = line2
				pair2 = line1				
	
			#--Tomamos los datos del diccionario y analizamos la lectura en la que se detectó el SL
			if sl[read_name][-3] == 1: #--SL detectado en el par 1
				col = pair1.split('\t')
				bitwise= int(col[1])
				bits = bin(bitwise)[2:]
				byte = "{0:0>12}".format(bits)
	
			else: #--SL detectado en el par 2
				col = pair2.split('\t')
				bitwise= int(col[1])
				bits = bin(bitwise)[2:]
				byte = "{0:0>12}".format(bits)
	
			#--The read do not align, not trim
			if int(byte[-3]) == 1:		
				del sl[col[0]]
				continue	
	
			#--The read align
			else:
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
						#--Escribimos las coordenadas del fragmento eliminido en un fichero gtf para extraer la secuencia genómica
						transcript = [col[2], 'CBMSO', 'transcript', str(x), str(y), '.', '+', '.']
						att = ['gene_id "', col[0], '"; ', 'transcript_id "',col[0], '"; '] 	
						gtf_temp.write('\t'.join(transcript) + '\t' + ''.join(att) + '\n')
	print('SAM analisis has finished...\n')	

	#--Extraemos las secuencias genómicas
	subprocess.call ('gffread -w ' + args.sample_id + '_temp.fa -g ' + genome + ' ' + args.sample_id + '_temp.gtf', shell = True)

	#--Reformateo del fichero temp.fa por si la secuencia ocupa de más de 1 línea (en el caso de que la secuencia eliminada sea mayor de 60 pb)
	refasta(args.sample_id + '_temp.fa')

	#--Analizamos las secuencias genómicas en busca de secuencias leader genómicas (mismas condiciones que en las lecturas)
	with open(args.sample_id + '_temp.fa', 'r') as seqs:
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
	with open (args.fqfile_r1, 'r') as fastq_infile1, open (args.fqfile_r2, 'r') as fastq_infile2, open(args.sample_id + '_SL_trimmed_1.fastq', 'w') as reads_outfile1, open(args.sample_id + '_SL_trimmed_2.fastq', 'w') as reads_outfile2, open(sl_seqs_file, 'w') as out_sl: #open('seqleader_ids', 'w') as leader_outfile
		count = 0

		while True:
			name1 = fastq_infile1.readline().rstrip('\n')
			seq1 = fastq_infile1.readline().rstrip('\n')
			coment1 = fastq_infile1.readline().rstrip('\n')
			qual1 = fastq_infile1.readline().rstrip('\n')

			name2 = fastq_infile2.readline().rstrip('\n')
			seq2 = fastq_infile2.readline().rstrip('\n')
			coment2 = fastq_infile2.readline().rstrip('\n')
			qual2 = fastq_infile2.readline().rstrip('\n')
	
			if not name1: #--Finish
				break
	
			name = name1.split(' ')[0][1:]
	
			if sl[name] and sl[name][-1]:
		
				count = count + 1
				new_name = name1.split(' ')[0] + 'L ' + ' '.join(name1.split(' ')[1:]) #--Añado una L al final del nombre de la lectura para identificarlas después
		
				if sl[name][-3] == 1: #--recortar el par 1		
					record1 = [new_name, seq1[sl[name][-2]:], coment1, qual1[sl[name][-2]:]]
					record2 = [new_name, seq2, coment2, qual2] 
					#leader_outfile.write(name + '\t' + str(sl[name][-3]) + '\n') #--Guardar sólo los ids y luego sacarlos del sam general
					reads_outfile1.write('\n'.join(record1) + '\n')
					reads_outfile2.write('\n'.join(record2) + '\n')
					out_sl.write('>' + new_name[1:] + '\n' + seq1[:sl[name][-2]] + '\n')
				else:
					record1 = [new_name, seq1, coment1, qual1]
					record2 = [new_name, seq2[sl[name][-2]:], coment2, qual2[sl[name][-2]:]] 
					#leader_outfile.write(name + '\t' + str(sl[name][-3]) + '\n') #--Guardar sólo los ids y luego sacarlos del sam general
					reads_outfile1.write('\n'.join(record1) + '\n')
					reads_outfile2.write('\n'.join(record2) + '\n')
					out_sl.write('>' + new_name[1:] + '\n' + seq2[:sl[name][-2]] + '\n')
			else:
				record1 = [name1, seq1, coment1, qual1]
				record2 = [name2, seq2, coment2, qual2]
				reads_outfile1.write('\n'.join(record1) + '\n')
				reads_outfile2.write('\n'.join(record2) + '\n')

	print(f'Valid leader hits: {count}\n')

	os.remove(args.sample_id + '_seqleader_temp_1.fastq')
	os.remove(args.sample_id + '_seqleader_temp_2.fastq')
	os.remove(args.sample_id + '_temp.sam')
	os.remove(args.sample_id + '_temp.gtf')
	os.remove(args.sample_id + '_temp.fa')


if __name__ == '__main__':
    main()


Valid leader hits: 925575

'''
