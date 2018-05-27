import random

import numpy
import pysam
from Crypto.Cipher import AES

count_reads = 0
count_candidates = 0
count_extensions_case_default = 0
count_extensions_case_reverse = 0
count_extensions_not_done = 0
count_contraction_case_default = 0
count_contraction_case_reverse = 0
count_contraction_not_done = 0
count_contraction_not_reverse_because_int = 0
count_contraction_not_reverse_because_string = 0
count_no_fingerprinting_because_equal=0


def byte_to_substring(byte):
    mask = 0b00000011
    values = ""
    #We read 4 nucleotides from a byte
    for i in range(0, 4):
        new_val = (byte & mask) >> (2 * i)
        if new_val == 1:
            values += 'A'
        elif new_val == 2:
            values += 'C'
        elif new_val == 3:
            values += 'T'
        else:
            values += 'G'
        mask <<= 2
    return values


def approximate_pdf(value_input):
    if value_input <= 166:
        encoded_value = 0
    elif value_input <= 189:
        encoded_value = 1
    elif value_input <= 202:
        encoded_value = 2
    elif value_input <= 208:
        encoded_value = 3
    else:
        encoded_value = value_input - 208 + 3
    return encoded_value


def random_nucletotide_gen():
    random_int = random.randrange(1, 5)
    if random_int == 1:
        return 'A'
    elif random_int == 2:
        return 'C'
    elif random_int == 3:
        return 'T'
    else:
        return 'G'


def perform_fingerprinting_start(read, read_s, reference, from_s, to_s, tentative_input_string, tentative_output_string):
    global count_extensions_case_default
    global count_extensions_case_reverse
    global count_extensions_not_done
    global count_no_fingerprinting_because_equal
    global count_contraction_case_default
    global count_contraction_case_reverse
    global count_contraction_not_done
    global count_contraction_not_reverse_because_int
    global count_contraction_not_reverse_because_string
    global count_no_fingerprinting_because_equal

    if from_s == to_s:
        count_no_fingerprinting_because_equal+=1
        return read, to_s

    if from_s < to_s:
        if from_s == read_s:
            result_transform = expand_s_start(read, from_s, to_s)
            result_s = to_s
            count_extensions_case_default += 1
        else:
            if to_s == read_s:
                result_transform = change_s_start(read, from_s)
                result_s = from_s
                count_extensions_case_reverse += 1
            else:
                result_transform = read
                result_s = read_s
                count_extensions_not_done += 1
    else:
        if from_s == read_s and read.query_sequence[:from_s] == tentative_input_string:
            result_transform = compact_s_start(read, from_s, to_s, tentative_output_string, reference)
            result_s = to_s
            count_contraction_case_default += 1
        else:
            if to_s == read_s and read.query_sequence[:to_s] == tentative_output_string:
                result_transform = compact_s_start(read, to_s, from_s, tentative_input_string, reference)
                result_s = from_s

                count_contraction_case_reverse += 1
            else:
                if to_s != read_s:
                    count_contraction_not_reverse_because_int += 1
                if read.query_sequence[to_s:from_s] != tentative_output_string:
                    count_contraction_not_reverse_because_string += 1
                result_transform = read
                result_s = read_s
                count_contraction_not_done += 1
    return result_transform, result_s


def undoing_fingerprinting(read, read_s, reference, from_s, to_s, tentative_input_string, tentative_output_string):
    if from_s < to_s:
        if read_s == to_s:
            undone_result = undo_expansion(read, reference, from_s, to_s)
            undone_result = change_s_start(read, from_s)
            undone_s = from_s
        else:
            if read_s == from_s:
                undone_result = change_s_start(read, to_s)
                undone_s = from_s
            else:
                undone_result = read
                undone_s = read_s
    else:
        if to_s == read_s and read.query_sequence[:to_s] == tentative_output_string:
            quality_backup = read.query_qualities
            new_query_sequence = tentative_input_string + read.query_sequence[from_s:]
            read.query_sequence = new_query_sequence
            undone_result = read
            change_s_start(read, from_s)
            undone_s = from_s
            read.query_qualities = quality_backup
        else:
            if from_s == read_s and read.query_sequence[:from_s] == tentative_input_string:
                quality_backup = read.query_qualities
                undone_query_sequence = tentative_output_string + reference[to_s:from_s]+read.query_sequence[from_s:]
                read.query_sequence = undone_query_sequence
                change_s_start(read, to_s)
                undone_s = to_s
                read.query_qualities = quality_backup
                undone_result = read
            else:
                undone_result = read
                undone_s = read_s

    return undone_result, undone_s


def change_s_start(read, new_s):
    cigar_tuples = read.cigartuples
    #Check if first operation is not a softclip
    if cigar_tuples[0][0] != 4:
	#next operation is shortened by the size of the softclip
        new_cigar = [(4, new_s)] + [(cigar_tuples[0][0], cigar_tuples[0][1] - new_s)] + cigar_tuples[1:]
    else:
        old_s = cigar_tuples[0][1]
	#we update the size of the softclip, and update the size of the next operation
        new_cigar = [(cigar_tuples[0][0], new_s)] + [
            (cigar_tuples[1][0], cigar_tuples[1][1] - (new_s - old_s))] + cigar_tuples[2:]
    filtered_cigar = [(type_cigar, count_cigar) for (type_cigar, count_cigar) in new_cigar if count_cigar != 0]
    read.cigartuples = filtered_cigar
    return read


def expand_s_start(read, old_s, new_s):
    quality_backup = read.query_qualities
    new_query_sequence = read.query_sequence[:old_s]
    #We obtain randomly the missing nucleotides
    for new_s_it in range(old_s, new_s):
        random_nucleotide = random_nucletotide_gen()
        while random_nucleotide == read.query_sequence[new_s_it]:
            random_nucleotide = random_nucletotide_gen()
        new_query_sequence += random_nucleotide
    new_query_sequence += read.query_sequence[new_s:]

    read = change_s_start(read, new_s)
    read.query_sequence = new_query_sequence
    read.query_qualities = quality_backup
    return read


def compact_s_start(read, old_s, new_s, substitution_string,reference):
    quality_backup = read.query_qualities
    new_query_sequence = substitution_string
    #We copy from the refence the missing nucleotides
    new_query_sequence += reference[new_s:old_s]
    new_query_sequence += read.query_sequence[max(old_s,new_s):]

    read = change_s_start(read, new_s)
    read.query_sequence = new_query_sequence
    read.query_qualities = quality_backup
    return read



def undo_expansion(read, reference, original_s, modified_s):
    quality_backup = read.query_qualities
    undone_query_sequence = read.query_sequence[:original_s]
    for new_match_it in range(original_s, modified_s):
        undone_query_sequence += reference[new_match_it]
    undone_query_sequence += read.query_sequence[modified_s:]
    read.query_sequence = undone_query_sequence
    read.query_qualities = quality_backup
    return read


class fingerprinter:
    max_I_in_a_row = 5
    max_D_in_a_row = 5
    max_N_in_a_row = 5
    max_IDN_in_read = 15

    class Entry:
        def __init__(self, cigar_tuples):
            self.is_candidate = True
            self.description = bytearray()
            self.S_start = 0
            self.S_end = 0

            if cigar_tuples[0][0] == 4:
                self.S_start = cigar_tuples[0][1]
            if cigar_tuples[-1][0] == 4:
                self.S_end = cigar_tuples[-1][1]

            type_i_in_read = 0
            type_d_in_read = 0
            type_n_in_read = 0
            type_s_in_read = 0
            current_position = 0
            for cigar_tuple in cigar_tuples:
                if cigar_tuple[0] != 0:
                    #Operation is Insert
                    if cigar_tuple[0] == 1:
                        type_i_in_read += cigar_tuple[1]
                        if cigar_tuple[1] > fingerprinter.max_I_in_a_row:
                            self.is_candidate = False
                            break
                        else:
                            self.description.append(current_position)
                            size_and_type_byte = cigar_tuple[1] << 4 | 0x00
                            self.description.append(size_and_type_byte)
                    #Operation is Delete
                    elif cigar_tuple[0] == 2:
                        type_d_in_read += cigar_tuple[1]
                        if cigar_tuple[1] > fingerprinter.max_D_in_a_row:
                            self.is_candidate = False
                            break
                        else:
                            self.description.append(current_position)
                            size_and_type_byte = cigar_tuple[1] << 4 | 0x01
                            self.description.append(size_and_type_byte)
                    #Operation is N
                    elif cigar_tuple[0] == 3:
                        type_n_in_read += cigar_tuple[1]
                        if cigar_tuple[1] > fingerprinter.max_N_in_a_row:
                            self.is_candidate = False
                            break
                        else:
                            self.description.append(current_position)
                            size_and_type_byte = cigar_tuple[1] << 4 | 0x02
                            self.description.append(size_and_type_byte)
                    elif cigar_tuple[0] != 4:
                        print("unknown cigar type: " + str(cigar_tuple[0]))
                current_position += cigar_tuple[1]
            if type_i_in_read + type_d_in_read + type_n_in_read > fingerprinter.max_IDN_in_read:
                self.is_candidate = False
            else:
                for i in range(len(self.description), 32):
                    self.description.append(0)

    @staticmethod
    def fingerprint(file_name):
        global count_reads
        global count_candidates
        global count_extensions_case_default
        global count_extensions_case_reverse
        global count_extensions_not_done
        global count_contraction_case_default
        global count_contraction_case_reverse
        global count_contraction_not_done
        global count_contraction_not_reverse_because_int
        global count_contraction_not_reverse_because_string
        global count_no_fingerprinting_because_equal

        AES_KEY = bytearray.fromhex('E961B2AC40AAC4CC36A8BF65BCA9177EA3CDE89190B7C35D1D259F452979D853')
        cipher = AES.new(bytes(AES_KEY), AES.MODE_ECB)

        file_name = file_name
        bamfile = pysam.AlignmentFile(file_name, "rb")
        fastafile = pysam.FastaFile("hs37d5.fa")

        fetched_reads = bamfile.fetch()

        count_fingerprint_candidate = 0

        for fetched_read in fetched_reads:
            count_reads += 1

            cigar_tuples = fetched_read.cigartuples
            entry = fingerprinter.Entry(cigar_tuples)
            if not entry.is_candidate:
                continue
            count_candidates += 1

            #We obtain the ciphertext associated to the generated description
            ciphertext = cipher.encrypt(bytes(entry.description))

            encoded_s_start = approximate_pdf(int(ciphertext[0]))
            encoded_new_s_start = approximate_pdf(int(ciphertext[1]))

            if encoded_s_start == encoded_new_s_start:
                continue

            encoded_s_end = approximate_pdf(int(ciphertext[2]))
            encoded_new_s_end = approximate_pdf(int(ciphertext[3]))

            if abs(encoded_new_s_start - encoded_s_start) > 10 or abs(encoded_new_s_end - encoded_s_end) > 10:
                #Arbitrary decision that the changes would be too significant
                continue

            cipher_position = 4
            tentative_original_string = ""
            tentative_new_string = ""

            #We read the encoded softclip operations from the ciphertext
            while len(tentative_original_string) < encoded_s_start:
                tentative_original_string = tentative_original_string + byte_to_substring(
                    ciphertext[cipher_position])
                cipher_position += 1
            tentative_original_string = tentative_original_string[:encoded_s_start]
            while len(tentative_new_string) < encoded_new_s_start:
                tentative_new_string = tentative_new_string + byte_to_substring(ciphertext[cipher_position])
                cipher_position += 1
            tentative_new_string = tentative_new_string[:encoded_new_s_start]


            count_modifiables = 0
            #We count as modifiable if the first two first operations are either an S or M
            if cigar_tuples[0][0] == 4:
                count_modifiables += cigar_tuples[0][1]
                if cigar_tuples[1][0] == 0:
                    count_modifiables += cigar_tuples[1][1]
            elif cigar_tuples[0][0] == 0:
                count_modifiables += cigar_tuples[0][1]

            if encoded_s_start < count_modifiables and encoded_new_s_start < count_modifiables:
                if (fetched_read.cigartuples[0][0] == 4):
                    offset_start_reference = fetched_read.cigartuples[0][1]
                else:
                    offset_start_reference = 0
                reference = fastafile.fetch(bamfile.get_reference_name(fetched_read.reference_id),
                                            fetched_read.reference_start - offset_start_reference,
                                            fetched_read.reference_start + fetched_read.infer_query_length())
                original_sequence = fetched_read.query_sequence
                original_cigar = fetched_read.cigartuples

                result_fingerprinting, fingerprinting_s = perform_fingerprinting_start(fetched_read, entry.S_start, reference,
                                                                                 encoded_s_start, encoded_new_s_start,
                                                                                 tentative_original_string,
                                                                                 tentative_new_string)
                copy_result_fingerprinting = result_fingerprinting.query_sequence
                result_undone, undone_s = undoing_fingerprinting(fetched_read, fingerprinting_s, reference, encoded_s_start,
                                                               encoded_new_s_start, tentative_original_string,
                                                               tentative_new_string)
                if result_undone.query_sequence != original_sequence or result_undone.cigartuples != original_cigar:
                    mismatch_indexes = [i for i in range(len(result_undone.query_sequence)) if result_undone.query_sequence[i] != original_sequence[i]]
                    problem_caused_by_reference = True
                    for mismatch_index in mismatch_indexes:
                        if result_undone.query_sequence[mismatch_index]!=reference[mismatch_index]:
                            problem_caused_by_reference= False
                            break
                    if not problem_caused_by_reference:
                        print("original: "+original_sequence)
                        print("fingerprint:"+copy_result_fingerprinting)
                        print("original: "+original_sequence)
                        print("undone:   "+result_undone.query_sequence)
                        print("reference:"+reference)
                        print("ori cigar:"+str(original_cigar))
                        print("undone ci:"+str(result_undone.cigartuples))
                        print("from: "+str(encoded_s_start)+" to: "+str(encoded_new_s_start)+ " encoded start: "+tentative_original_string+" encoded output: "+tentative_new_string)
                        print("ALERT")


            count_fingerprint_candidate += 1
        print("count_reads: " + str(count_reads))
        print("count_candidates: " + str(count_candidates))
        print("count_fingerprint_candidate: " + str(count_fingerprint_candidate))
        print("count_extensions_case_default: " + str(count_extensions_case_default))
        print("count_extensions_case_reverse: " + str(count_extensions_case_reverse))
        print("count_extensions_not_done: " + str(count_extensions_not_done))
        print("count_contraction_case_default: " + str(count_contraction_case_default))
        print("count_contraction_case_reverse: " + str(count_contraction_case_reverse))
        print("count_contraction_not_done: " + str(count_contraction_not_done))
        print("count_contraction_not_reverse_because_int: " + str(count_contraction_not_reverse_because_int))
        print("count_contraction_not_reverse_because_string: " + str(count_contraction_not_reverse_because_string))
        print("count_no_fingerprinting_because_equal: "+str(count_no_fingerprinting_because_equal))
