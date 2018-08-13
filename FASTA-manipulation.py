path = 'protein_fasta1.fasta'

def read_fasta_file(filepath):
    fasta_file = open(filepath,'r')
    #fasta_file = fasta_file.read().replace('\n', '')
    return fasta_file

prepped_fasta = read_fasta_file(path)

def parse_fasta_record(fasta):
    fasta_dict = {'Accession':'','Sequence':'','Metadata':''}
    fasta_read = fasta.read()
    #print fasta_read
    fasta_lines = fasta_read.splitlines()
    #print fasta_lines
    fasta_stripped = fasta_read.replace('\n', '')
    #print fasta_stripped
    split_fasta = fasta_stripped.split(">")[1:]
    #print split_fasta
    for line in fasta_lines:
        if line.startswith(">"):
            meta = line[17:]
            #print meta
            fasta_dict['Metadata'] += meta
        if line[0].isalpha():
            seq = line[:]
            #print seq
            fasta_dict['Sequence'] += seq
    for n in split_fasta:
        acc = n[1:14].replace('\n','')
        fasta_dict['Accession'] += acc
    return fasta_dict
       
#print parse_fasta_record(prepped_fasta)
parsed_dict = parse_fasta_record(prepped_fasta)
print parsed_dict
print "****"

#Find seq_type
if "protein" in parsed_dict['Metadata']:
    seq_type = "Protein"
elif "DNA" in parsed_dict['Metadata']:
    seq_type = "DNA"
elif "RNA" in parsed_dict['Metadata']:
    seq_type = "RNA"
else:
    print "Not found"
seq_type = seq_type
#print type(seq_type)

# from itertools import groupby
import random
class FASTA_File:
    def __init__(self, filepath, seq):
        self.filepath = filepath
        self.seq = seq
    def fasta_dictionary(self, filepath, seq):
        fasta_dict = {'Sequence':'', 'Metadata':''}
        fasta_dict['Sequence'] += seq
        fasta_dict['Metadata'] += parsed_dict['Metadata']
        half_dict = fasta_dict
        key = parsed_dict['Accession']
        self.fasta_dictionary = {key:fasta_dict}
        return self.fasta_dictionary
    def accession_random_sample(self, filepath, seq):
        iso_acc = parsed_dict['Accession']
        acc_split = iso_acc.split(" ")[:-1]
        sample_number = random.randint(1,len(acc_split))
        ran_sam = random.sample(population=acc_split, k=sample_number)
        return ran_sam
    def __and__(self, other):
        dict_acc = parsed_dict['Accession']
        split_it = dict_acc.split(" ")[:-1]
        rand1 = self.accession_random_sample(path, split_it)
        rand2 = self.accession_random_sample(path, split_it)
        #print rand1
        #print rand2
        inter_list = []
        for int in rand1:
            if int in rand2:
                return int
    def __or__(self, other):
        dict_acc = parsed_dict['Accession']
        split_it = dict_acc.split(" ")[:-1]
        rand3 = self.accession_random_sample(path, split_it)
        int1 = self.__and__(split_it)
        #print rand3
        #print int1
        union = rand3 + [int1]
        return union
    def subset_write_out(self, accessions):
        self.accessions = accessions
        #print self.accessions
        #print parsed_dict
        if parsed_dict['Accession'] == self.accessions:
            with open(path,'r') as f:
                return f.read()
    def in_metadata(self, search):
        self.search = search#.lower() #for non-case sensitive search
        dict_meta = parsed_dict['Metadata']
        dict_acc = parsed_dict['Accession']
        acc_list = dict_acc.split(" ")[:-1]
        if self.search in dict_meta:
            #use dict_meta.lower() for non-case sensitive search
            return acc_list
        else:
            return "Not found"
        
class Sequence:
    def __init__(self, filepath):
        self.filepath = filepath
    def randomized_sequence(self, filepath):
        iso_seq = parsed_dict['Sequence']
        lst_seq = list(iso_seq)
        random.shuffle(lst_seq)
        ran_seq = ''.join(lst_seq)
        pick_seq = ''.join(random.choice(ran_seq) for x in range(random.randint(1,len(iso_seq)))) 
        #print len(ran_seq)
        #print len(pick_seq)
        return pick_seq
class Nucleotide:
    def __init__(self, subsequence):
        self.subsequence = subsequence
    def mask(self, subsequence):
        #print self.subsequence
        lower_sub = self.subsequence.lower()
        dict_seq = parsed_dict['Sequence']
        if self.subsequence in dict_seq:
            return dict_seq.replace(self.subsequence,lower_sub)

#Print statements        
FASTA_mani = FASTA_File(path, seq_type)
print FASTA_mani.fasta_dictionary(path, seq_type)
print "****"
print FASTA_mani.accession_random_sample(path, seq_type)
print "****"
print FASTA_mani.__and__(path)
print "****"
print FASTA_mani.__or__(path)
print "****"
print FASTA_mani.subset_write_out(parsed_dict['Accession'])
print "****"
print FASTA_mani.in_metadata("length=1434")
# test "PlaSMODium falCiparum", "length=1434", "pyRUVate kinAse"
print "********"
Seq_mani = Sequence(path)
print Seq_mani.randomized_sequence(path)
print "********"
Nuc_mani = Nucleotide("KKIYNK")
print Nuc_mani.mask("KKIYNK")
# test "NNN", "KKIYNK", "WKKGEDLIFLEEERDKELNLQIEK"


#Unit test
if __name__=="__main__":
    import unittest
    
    def test_read_fasta_file():
        read_fasta = read_fasta_file(path)
        answer = prepped_fasta
        assertEqual(read_fasta.read_fasta_file(),answer)
        
    def test_parse_fasta_record():
        parse_fasta = parse_fasta_record(prepped_fasta)
        answer = parsed_dict
        assertEqual(parse_fasta.parse_fasta_record(),answer)
    
    class Test_FASTA_File_Objects(unittest.TestCase):
        def test_fasta_dictionary(self):
            FASTA_mani = FASTA_File(path, seq_type)
            fasta_dict = FASTA_mani.fasta_dictionary(path, seq_type)
            dict_type = type(fasta_dict)
            answer = type({})
            if dict_type == answer:
                return True
        def test__and__(self):
            FASTA_mani = FASTA_File(path, seq_type)
            inter = str(FASTA_mani.__and__(path))
            if inter in parsed_dict['Accession']:
                return True
        def test__or__(self):
            FASTA_mani = FASTA_File(path, seq_type)
            uni = FASTA_mani.__or__(path)
            for u in uni:
                if u in parsed_dict['Accession']:
                    return True
        def in_metadata(self):
            FASTA_mani = FASTA_File(path, seq_type)
            meta = FASTA_mani.in_metadata("length=1434")
            answer = parsed_dict['Accession']
            self.assertEqual(meta.in_metadata(),answer)
            
    class Test_Nucleotide_Objects(unittest.TestCase):
        def test_mask(self):
            mask_result = Nuc_mani.mask("KKIYNK")
            answer = parsed_dict['Sequence']
            if mask_result != answer:
                return True

    unittest.main()