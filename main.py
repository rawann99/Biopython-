import getopt
import sys
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio.pairwise2 import format_alignment, align
from Bio import pairwise2
from Bio.SeqUtils import *
import fileinput
import glob
from Bio import SeqIO
se='TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG'
result_handle = NCBIWWW.qblast("blastn", "nt", se)
blast_record = NCBIXML.read(result_handle)
for alignment in blast_record.descriptions:
    print("Tittle", alignment.title)



def gc(seq):
    if is_valid(seq, "dna"):
        sequence = Seq(seq)
        print(GC(sequence))
    else:
        print("ENTER VALID SEQ")


def transcrib(seq):
    if is_valid(seq, "dna"):
        sequence = Seq(seq)
        print(sequence.transcribe())
    else:
        print("ENTER VALID SEQ")


def reverse_complement(seq):
    if is_valid(seq, "dna"):
        sequence = Seq(seq)
        print(sequence.reverse_complement())
    else:
        print("ENTER VALID SEQ")


def calc_nbases(seq):
    ncounter = seq.count("n")
    ncounter += seq.count("N")
    print(ncounter)


def is_valid(seq, type):
    type = type.lower()
    valid = True
    seq = seq.upper()
    sequences = {'protein': 'SDVIHRYKUUPAKSHGWYVCJRSRFTWMVWWRFRSCRA', 'dna': 'AGCT', 'rna': 'AGCU'}
    if type not in sequences.keys():
        print("wrong type")
        return False
    for i in seq:
        if i not in sequences[type]:
            valid = False
    return valid


def filter_nbases(seq):
    correct = ''
    seq = seq.upper()
    for i in seq:
        if i not in "ATCG":
            continue
        else:
            correct += i
    print(correct)


def seq_alignment(seq1, seq2, file=''):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    seqq1 = Seq(seq1)
    seqq2 = Seq(seq2)
    alignments = pairwise2.align.globalxx(seqq1, seqq2)
    if file == '':
        for alignment in alignments:
            print(format_alignment(*alignment))
    else:
        f = open(file, "w")
        for alignment in alignments:
            print(format_alignment(*alignment))
            f.write(format_alignment(*alignment))
        f.close()


def seq_alignment_files(file1, file2, file=''):
    type_file_1 = file1.split('.')[1]
    type_file_2 = file2.split('.')[1]
    seq1 = SeqIO.parse(file1, type_file_1)
    seq2 = SeqIO.parse(file2, type_file_2)
    for i, j in zip(seq1, seq2):
        alignments = pairwise2.align.globalxx(i.seq, j.seq)
        if file == '':
            for alignment in alignments:
                print(format_alignment(*alignment))
        else:
            f = open(file, "a")
            for alignment in alignments:
                print(format_alignment(*alignment))
                f.write(format_alignment(*alignment))
            f.close()


def online_alignment(seq, file=''):
    if file == '':
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)
        blast_record = NCBIXML.read(result_handle)
        print("Len:", len(blast_record.alignments))
        for alignment in blast_record.descriptions:
            print("Tittle", alignment.title)
            print("Score", alignment.score)
            print("e", alignment.e)
            print("num_alignments", alignment.num_alignments)
        for alignment in blast_record.alignments:
            print("Tittle", alignment.title)
            print("Length", alignment.length)
            for hsp in alignment.hsps:
                print("HSP expect", hsp.expect)
                print("HSP score", hsp.score)
                print("HSP bits", hsp.bits)
                print("HSP query", hsp.query)
                print("HSP match", hsp.match)
                print("HSP sbjct", hsp.sbjct)
                print("HSP gaps", hsp.gaps)
                print("HSP.num_alignments", hsp.num_alignments)
                print("HSP.positives", hsp.positives)
    else:
        f = open(file, "a")
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)
        blast_record = NCBIXML.read(result_handle)
        f.write(str(len(blast_record.alignments)))
        for alignment in blast_record.descriptions:
            print("Tittle", alignment.title)
            f.write(alignment.title)
            print("Score", alignment.score)
            f.write(str(alignment.score))
            print("e", alignment.e)
            f.write(str(alignment.e))
            print("num_alignments", alignment.num_alignments)
            f.write(str(alignment.num_alignments))
        for alignment in blast_record.alignments:
            print("Tittle", alignment.title)
            f.write(alignment.title)
            print("Length", alignment.length)
            f.write(str(alignment.length))
            for hsp in alignment.hsps:
                print("HSP expect", hsp.expect)
                f.write(str(hsp.expect))
                print("HSP score", hsp.score)
                f.write(str(hsp.score))
                print("HSP bits", hsp.bits)
                f.write(str(hsp.bits))
                print("HSP query", hsp.query)
                f.write(hsp.query)
                print("HSP match", hsp.match)
                f.write(str(hsp.match))
                print("HSP sbjct", hsp.sbjct)
                f.write(hsp.sbjct)
                print("HSP gaps", hsp.gaps)
                f.write(str(hsp.gaps))
                print("HSP.num_alignments", hsp.num_alignments)
                f.write(str(hsp.num_alignments))
                print("HSP.positives", hsp.positives)
                f.write(str(hsp.positives))

        f.close()


def merge_fasta(*files, file=''):
    for i in files:
        if ".fasta" not in i:
            print("YOU MUST ENTER FASTA FILES!!")
            return False
    if file == '':
        for i in files:
            input_lines = fileinput.input(i)
            for j in input_lines:
                print(str(j))

    else:
        with open(file, 'a') as file:
            for i in files:
                input_lines = fileinput.input(i)
                file.writelines(input_lines)
            file.close()


def convert_to_fasta(file):
    if ".gb" in file:
        with open(file) as input_handle, open("OUT.fasta", "w") as output_handle:
            sequences = SeqIO.parse(input_handle, "genbank")
            SeqIO.write(sequences, output_handle, "fasta")
    else:
        print("YOU MUST ENTER FILE.gb extension")
merge_fasta('fasta1.fasta','fasta2.fasta','fasta3.fasta','dar4.fasta')

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 'o:')
    # opts contain all options with args || args contain all args that have no optin
    print("opts:",opts)
    print(type(opts[0]))
    print(sys.argv[0])
    print("args:",args)
    if args[0] == 'gc':
        gc(args[1])

    elif args[0] == 'transcrib':
        transcrib(args[1])

    elif args[0] == 'reverse_complement':
        reverse_complement(args[1])

    elif args[0] == 'calc_nbases':
        calc_nbases(args[1])

    elif args[0] == 'is_valid':
        print(is_valid(args[1], args[2]))

    elif args[0] == 'filter_nbases':
        filter_nbases(args[1])

    elif args[0] == 'seq_alignment':
        if len(args[1:]) == 2:
            if len(opts) != 0:
                if opts[0][0] == '-o':
                    seq_alignment(args[1], args[2], opts[0][1])
            else:
                seq_alignment(args[1], args[2])
        else:
            print("YOU MUST ENTER 2 ARGUMENT ONLY!!")

    elif args[0] == 'seq_alignment_files':
        if len(args[1:]) == 2:
            if len(opts) != 0:
                if opts[0][0] == '-o':
                    seq_alignment_files(args[1], args[2], opts[0][1])
            else:
                seq_alignment_files(args[1], args[2])
        else:
            print("YOU MUST ENTER 2 ARGUMENT ONLY!!")

    elif args[0] == 'online_alignment':
        if len(args[1:]) == 1:
            if len(opts) != 0:
                if opts[0][0] == '-o':
                    online_alignment(args[1], opts[0][1])
            else:
                online_alignment(args[1])
        else:
            print("YOU MUST ENTER 2 ARGUMENT ONLY!!")

    elif args[0] == "merge_fasta":
        if len(args[1:]) < 2:
            print("YOU MUST ENTER AT LEAST 2 ARGUMENT!!")
        else:
            if len(opts) != 0:
                if opts[0][0] == '-o':
                    tup = tuple(args[1:])
                    merge_fasta(tup, file=opts[0][1])

            else:
                tup = tuple(args[1:])
                merge_fasta(tup)

    elif args[0] == "convert_to_fasta":
        if len(args[1:]) == 1:

            convert_to_fasta(args[1])
        else:
            print("YOU MUST ENTER ONLY ONE FILE!!")

except getopt.GetoptError as err:
    print(err)