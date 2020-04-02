from __future__ import print_function
'''
Created on Dec 31, 2011

@author: dominic
'''
import sys, os, subprocess
sys.path.append("../")


class FormatDB():
    def __init__(self, _in, fdb_path, input_type='fasta', dbtype = 'nucl', title = ""):
        self.fdb_path = fdb_path
        self.parse_seqids = False
        self.ud = os.path.dirname(fdb_path)
        #pathlib.Path(self.ud).mkdir(parents=True, exist_ok=True)
        if not(os.path.exists(self.ud)):os.mkdir(self.ud)
        self.d = {'in':_in, 
                  'out':self.fdb_path, 
                  'input_type':input_type,
                  'dbtype':dbtype,
                  }
        self.title = title
        self.cmd = ""
        self.seq_len_dict ={}
        self.goup_sz_dict = {}
        self.n = 0
        if not(_in):self.compute_stats()
    def make_cmd(self):
        args = ' '.join(["-%s %s"%(k, self.d[k]) for k in self.d.keys()])
        self.cmd = self.cmd + args
    def RunFromFDB(self, parse_seqids=False):
        self.parse_seqids = parse_seqids
        self.make_cmd()
        print ('cmdstr:', "cd %s;makeblastdb  %s -parse_seqids;"%(self.ud, self.cmd))
        if parse_seqids:
            os.system("cd %s;makeblastdb  %s -parse_seqids;"%(self.ud, self.cmd))
        else:
            os.system("cd %s;makeblastdb  %s;"%(self.ud, self.cmd))
        self.compute_stats()
        return os.path.join(self.ud, self.d['out'])
    def AlreadyExists(self):
        test_path = self.d['out']+'.nsq'
        v = False
        if os.path.exists(test_path): 
            v = True
        print ('FormatDB::AlreadyExists: %s  exists: %s'%(test_path, str(v)))
        return v
    def compute_stats(self):
        self.seq_len_dict = self.MakeSequenceLengthLst()
        self.goup_sz_dict = self.MakeGroupSzLst()
        self.n = len(self.seq_len_dict)
    def ListBlastDB(self, Verbose = None):
        p = os.path.join(os.path.dirname(self.fdb_path), 'SLST.TXT')
        if self.parse_seqids:
            #using this format: >lcl|ACCESSION_ID     . . . and nothing else!!
            #https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastdbcmd_application_opti/
            cmdstr = "blastdbcmd -db %s -entry all "%(self.fdb_path) + '-outfmt "%a %l"' +  '|tee %s'%p
        else:
            cmdstr = "blastdbcmd -db %s -entry all "%(self.fdb_path) + '-outfmt "%t %l"' +  '|tee %s'%p
        print ('cmdstr:', cmdstr)
        os.system(cmdstr)
        with open(p, 'r') as f:
            result = f.read()
        lst = [r.replace('"', '').strip() for r in result.split('\n')]
        olst = []
        for rec in lst:
            if not(rec):continue
            olst.append(rec)
        print ('ListBlastDB: returning %i recs'%len(olst))
        return [f.split(' ') for f in olst]
    def MakeSequenceLengthLst(self, Verbose = None):
        lenLst = self.ListBlastDB()
        #lenLst = [f.split(' ') for f in faLst]
        seq_len_dict = {}
        for g in lenLst:
            if Verbose:print (g)
            seq_len_dict[g[0]] = int(g[-1])
        return seq_len_dict
    def MakeGroupSzLst(self):
        faLst = self.ListBlastDB()
        groupLst = [f[0].split('|')[0] for f in faLst]
        group_sz_dict = {}
        for g in groupLst:
            if not g in group_sz_dict: group_sz_dict[g] = 0.0
            group_sz_dict[g] += 1.0
        return group_sz_dict
    def GetSeqInterval(self, seq_name):
        #In order for this to work you need to use parse_seqids and use the following formats:
        #https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.T5
        #Options: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastdbcmd_application_opti/
        
        p = os.path.join(os.path.dirname(self.fdb_path), 'SLST.TXT')
        cmdstr = "blastdbcmd -db %s -entry %s "%(self.fdb_path, seq_name) + '-outfmt "%o %t %l"' +  '|tee %s'%p
        print ('cmdstr:', cmdstr)
        os.system(cmdstr)
        with open(p, 'r') as f:
            result = f.read()
        return result
    def __str__(self):
        return "FormatDB %s %i seqs"%(self.fdb_path, self.n)
        

class BlastRunObjOld():
    #TODO FIX THIS
    #Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output
    def __init__(self, query, _out, db, dust="no", outfmt="7", num_threads="8", word_size="16", num_descriptions="10000", rank=None):
        self.d = {'query':query, 
                  'out':_out, 
                  'db':db,
                  'dust':dust,
                  'outfmt':outfmt,
                  'num_threads':num_threads,
                  'word_size':word_size,
                  'max_target_seqs':num_descriptions
                  }
        self.cmd = ""
        self.rank = rank
        self.make_cmd()
    def make_cmd(self):
        args = ' '.join(["-%s %s"%(k, str(self.d[k])) for k in self.d.keys()])
        self.cmd = "blastn  %s;"%(args)
    def RunBlast(self):
        self.make_cmd()
        if self.rank:
            print ("[%i]"%self.rank,)
        print (self.cmd)
        os.system(self.cmd)
        return self.d['out']
    
class BlastRunObj():
    """Runs BLAST"""
    #For remote you need to install the latest version of blast ncbi-blast-2.8.1+-x64-linux.tar.gz
    #ftp://ftp.ncbi.nlm.nih.gov/pub/factsheets/HowTo_BLASTGuide.pdf <-- Guide
    """blastn [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-negative_seqidlist filename]
    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
    [-negative_taxidlist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-perc_identity float_value] [-qcov_hsp_perc float_value]
    [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value] [-penalty penalty]
    [-reward reward] [-no_greedy] [-min_raw_gapped_score int_value]
    [-template_type type] [-template_length int_value] [-dust DUST_options]
    [-filtering_db filtering_database]
    [-window_masker_taxid window_masker_taxid]
    [-window_masker_db window_masker_db] [-soft_masking soft_masking]
    [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
    [-best_hit_score_edge float_value] [-subject_besthit]
    [-window_size int_value] [-off_diagonal_range int_value]
    [-use_index boolean] [-index_name string] [-lcase_masking]
    [-query_loc range] [-strand strand] [-parse_deflines] [-outfmt format]
    [-show_gis] [-num_descriptions int_value] [-num_alignments int_value]
    [-line_length line_length] [-html] [-max_target_seqs num_sequences]
    [-num_threads int_value] [-remote] [-version]
    
    DESCRIPTION
   Nucleotide-Nucleotide BLAST 2.8.1+
        OPTIONAL ARGUMENTS
         -h
           Print USAGE and DESCRIPTION;  ignore all other parameters
         -help
           Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
         -version
           Print version number;  ignore other arguments

         *** Input query options
         -query <File_In>
           Input file name
           Default = `-'
         -query_loc <String>
           Location on the query sequence in 1-based offsets (Format: start-stop)
         -strand <String, `both', `minus', `plus'>
           Query strand(s) to search against database/subject
           Default = `both'

         *** General search options
         -task <String, Permissible values: 'blastn' 'blastn-short' 'dc-megablast'
                        'megablast' 'rmblastn' >
           Task to execute
           Default = `megablast'
         -db <String>
           BLAST database name
            * Incompatible with:  subject, subject_loc
         -out <File_Out>
           Output file name
           Default = `-'
         -evalue <Real>
           Expectation value (E) threshold for saving hits
           Default = `10'
         -word_size <Integer, >=4>
           Word size for wordfinder algorithm (length of best perfect match)
         -gapopen <Integer>
           Cost to open a gap
         -gapextend <Integer>
           Cost to extend a gap
         -penalty <Integer, <=0>
           Penalty for a nucleotide mismatch
         -reward <Integer, >=0>
           Reward for a nucleotide match
         -use_index <Boolean>
           Use MegaBLAST database index
           Default = `false'
         -index_name <String>
           MegaBLAST database index name (deprecated; use only for old style indices)

         *** BLAST-2-Sequences options
         -subject <File_In>
           Subject sequence(s) to search
            * Incompatible with:  db, gilist, seqidlist, negative_gilist,
           negative_seqidlist, taxids, taxidlist, negative_taxids, negative_taxidlist,
           db_soft_mask, db_hard_mask
         -subject_loc <String>
           Location on the subject sequence in 1-based offsets (Format: start-stop)
            * Incompatible with:  db, gilist, seqidlist, negative_gilist,
           negative_seqidlist, taxids, taxidlist, negative_taxids, negative_taxidlist,
           db_soft_mask, db_hard_mask, remote

         *** Formatting options
         -outfmt <String>
           alignment view options:
             0 = Pairwise,
             1 = Query-anchored showing identities,
             2 = Query-anchored no identities,
             3 = Flat query-anchored showing identities,
             4 = Flat query-anchored no identities,
             5 = BLAST XML,
             6 = Tabular,
             7 = Tabular with comment lines,
             8 = Seqalign (Text ASN.1),
             9 = Seqalign (Binary ASN.1),
            10 = Comma-separated values,
            11 = BLAST archive (ASN.1),
            12 = Seqalign (JSON),
            13 = Multiple-file BLAST JSON,
            14 = Multiple-file BLAST XML2,
            15 = Single-file BLAST JSON,
            16 = Single-file BLAST XML2,
            17 = Sequence Alignment/Map (SAM),
            18 = Organism Report

           Options 6, 7, 10 and 17 can be additionally configured to produce
           a custom format specified by space delimited format specifiers.
           The supported format specifiers for options 6, 7 and 10 are:
                    qseqid means Query Seq-id
                       qgi means Query GI
                      qacc means Query accesion
                   qaccver means Query accesion.version
                      qlen means Query sequence length
                    sseqid means Subject Seq-id
                 sallseqid means All subject Seq-id(s), separated by a ';'
                       sgi means Subject GI
                    sallgi means All subject GIs
                      sacc means Subject accession
                   saccver means Subject accession.version
                   sallacc means All subject accessions
                      slen means Subject sequence length
                    qstart means Start of alignment in query
                      qend means End of alignment in query
                    sstart means Start of alignment in subject
                      send means End of alignment in subject
                      qseq means Aligned part of query sequence
                      sseq means Aligned part of subject sequence
                    evalue means Expect value
                  bitscore means Bit score
                     score means Raw score
                    length means Alignment length
                    pident means Percentage of identical matches
                    nident means Number of identical matches
                  mismatch means Number of mismatches
                  positive means Number of positive-scoring matches
                   gapopen means Number of gap openings
                      gaps means Total number of gaps
                      ppos means Percentage of positive-scoring matches
                    frames means Query and subject frames separated by a '/'
                    qframe means Query frame
                    sframe means Subject frame
                      btop means Blast traceback operations (BTOP)
                    staxid means Subject Taxonomy ID
                  ssciname means Subject Scientific Name
                  scomname means Subject Common Name
                sblastname means Subject Blast Name
                 sskingdom means Subject Super Kingdom
                   staxids means unique Subject Taxonomy ID(s), separated by a ';'
                                 (in numerical order)
                 sscinames means unique Subject Scientific Name(s), separated by a ';'
                 scomnames means unique Subject Common Name(s), separated by a ';'
                sblastnames means unique Subject Blast Name(s), separated by a ';'
                                 (in alphabetical order)
                sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
                                 (in alphabetical order)
                    stitle means Subject Title
                salltitles means All Subject Title(s), separated by a '<>'
                   sstrand means Subject Strand
                     qcovs means Query Coverage Per Subject
                   qcovhsp means Query Coverage Per HSP
                    qcovus means Query Coverage Per Unique Subject (blastn only)
           When not provided, the default value is:
           'qaccver saccver pident length mismatch gapopen qstart qend sstart send
           evalue bitscore', which is equivalent to the keyword 'std'
           The supported format specifier for option 17 is:
                        SQ means Include Sequence Data
                        SR means Subject as Reference Seq
           Default = `0'
         -show_gis
           Show NCBI GIs in deflines?
         -num_descriptions <Integer, >=0>
           Number of database sequences to show one-line descriptions for
           Not applicable for outfmt > 4
           Default = `500'
            * Incompatible with:  max_target_seqs
         -num_alignments <Integer, >=0>
           Number of database sequences to show alignments for
           Default = `250'
            * Incompatible with:  max_target_seqs
         -line_length <Integer, >=1>
           Line length for formatting alignments
           Not applicable for outfmt > 4
           Default = `60'
         -html
           Produce HTML output?

         *** Query filtering options
         -dust <String>
           Filter query sequence with DUST (Format: 'yes', 'level window linker', or
           'no' to disable)
           Default = `20 64 1'
         -filtering_db <String>
           BLAST database containing filtering elements (i.e.: repeats)
         -window_masker_taxid <Integer>
           Enable WindowMasker filtering using a Taxonomic ID
         -window_masker_db <String>
           Enable WindowMasker filtering using this repeats database.
         -soft_masking <Boolean>
           Apply filtering locations as soft masks
           Default = `true'
         -lcase_masking
           Use lower case filtering in query and subject sequence(s)?

         *** Restrict search or results
         -gilist <String>
           Restrict search of database to list of GIs
            * Incompatible with:  seqidlist, taxids, taxidlist, negative_gilist,
           negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
           subject_loc
         -seqidlist <String>
           Restrict search of database to list of SeqIDs
            * Incompatible with:  gilist, taxids, taxidlist, negative_gilist,
           negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
           subject_loc
         -negative_gilist <String>
           Restrict search of database to everything except the specified GIs
            * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
           negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
           subject_loc
         -negative_seqidlist <String>
           Restrict search of database to everything except the specified SeqIDs
            * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
           negative_gilist, negative_taxids, negative_taxidlist, remote, subject,
           subject_loc
         -taxids <String>
           Restrict search of database to include only the specified taxonomy IDs
           (multiple IDs delimited by ',')
            * Incompatible with:  gilist, seqidlist, taxidlist, negative_gilist,
           negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
           subject_loc
         -negative_taxids <String>
           Restrict search of database to everything except the specified taxonomy IDs
           (multiple IDs delimited by ',')
            * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
           negative_gilist, negative_seqidlist, negative_taxidlist, remote, subject,
           subject_loc
         -taxidlist <String>
           Restrict search of database to include only the specified taxonomy IDs
            * Incompatible with:  gilist, seqidlist, taxids, negative_gilist,
           negative_seqidlist, negative_taxids, negative_taxidlist, remote, subject,
           subject_loc
         -negative_taxidlist <String>
           Restrict search of database to everything except the specified taxonomy IDs
            * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
           negative_gilist, negative_seqidlist, negative_taxids, remote, subject,
           subject_loc
         -entrez_query <String>
           Restrict search with the given Entrez query
            * Requires:  remote
         -db_soft_mask <String>
           Filtering algorithm ID to apply to the BLAST database as soft masking
            * Incompatible with:  db_hard_mask, subject, subject_loc
         -db_hard_mask <String>
           Filtering algorithm ID to apply to the BLAST database as hard masking
            * Incompatible with:  db_soft_mask, subject, subject_loc
         -perc_identity <Real, 0..100>
           Percent identity
         -qcov_hsp_perc <Real, 0..100>
           Percent query coverage per hsp
         -max_hsps <Integer, >=1>
           Set maximum number of HSPs per subject sequence to save for each query
         -culling_limit <Integer, >=0>
           If the query range of a hit is enveloped by that of at least this many
           higher-scoring hits, delete the hit
            * Incompatible with:  best_hit_overhang, best_hit_score_edge
         -best_hit_overhang <Real, (>0 and <0.5)>
           Best Hit algorithm overhang value (recommended value: 0.1)
            * Incompatible with:  culling_limit
         -best_hit_score_edge <Real, (>0 and <0.5)>
           Best Hit algorithm score edge value (recommended value: 0.1)
            * Incompatible with:  culling_limit
         -subject_besthit
           Turn on best hit per subject sequence
         -max_target_seqs <Integer, >=1>
           Maximum number of aligned sequences to keep
           (value of 5 or more is recommended)
           Default = `500'
            * Incompatible with:  num_descriptions, num_alignments

         *** Discontiguous MegaBLAST options
         -template_type <String, `coding', `coding_and_optimal', `optimal'>
           Discontiguous MegaBLAST template type
            * Requires:  template_length
         -template_length <Integer, Permissible values: '16' '18' '21' >
           Discontiguous MegaBLAST template length
            * Requires:  template_type

         *** Statistical options
         -dbsize <Int8>
           Effective length of the database
         -searchsp <Int8, >=0>
           Effective length of the search space

         *** Search strategy options
         -import_search_strategy <File_In>
           Search strategy to use
            * Incompatible with:  export_search_strategy
         -export_search_strategy <File_Out>
           File name to record the search strategy used
            * Incompatible with:  import_search_strategy

         *** Extension options
         -xdrop_ungap <Real>
           X-dropoff value (in bits) for ungapped extensions
         -xdrop_gap <Real>
           X-dropoff value (in bits) for preliminary gapped extensions
         -xdrop_gap_final <Real>
           X-dropoff value (in bits) for final gapped alignment
         -no_greedy
           Use non-greedy dynamic programming extension
         -min_raw_gapped_score <Integer>
           Minimum raw gapped score to keep an alignment in the preliminary gapped and
           traceback stages
         -ungapped
           Perform ungapped alignment only?
         -window_size <Integer, >=0>
           Multiple hits window size, use 0 to specify 1-hit algorithm
         -off_diagonal_range <Integer, >=0>
           Number of off-diagonals to search for the 2nd hit, use 0 to turn off
           Default = `0'

         *** Miscellaneous options
         -parse_deflines
           Should the query and subject defline(s) be parsed?
         -num_threads <Integer, >=1>
           Number of threads (CPUs) to use in the BLAST search
           Default = `1'
            * Incompatible with:  remote
         -remote
           Execute search remotely?
            * Incompatible with:  gilist, seqidlist, taxids, taxidlist,
           negative_gilist, negative_seqidlist, negative_taxids, negative_taxidlist,
           subject_loc, num_threads

    """
    def __init__(self, query, _out, db, blast_params={}):
        self.d = {'query':query, 
                  'out':_out, 
                  'db':db,
                  'dust':"no",
                  'outfmt':7,
                  'num_threads':8,
                  'word_size':17,
                  'max_target_seqs': 100000
                  }
        self.d.update(blast_params)
        ##Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output
        if db=="nr":
            self.d["remote"] = ""
            del self.d['num_threads']
        self.cmd = ""
        self.make_cmd(db=db)
    def make_cmd(self, db=None):
        args = ' '.join(["-%s %s"%(k, str(self.d[k])) for k in self.d.keys()])
        self.cmd = "blastn  %s;"%(args)
    def RunBlast(self, Verbose=False):
        self.make_cmd()
        if Verbose:print(self.cmd)
        ret = os.system(self.cmd)
        if Verbose:print("return from blast_run:", ret)
        return self.d['out']
    def Blast2Seqs(self, subj_path, out_path, params={}, Verbose=False):
        self.d.update(params)
        if 'db' in self.d: del self.d['db']
        self.d['subject'] = subj_path
        del self.d['num_threads']
        self.d['out'] = out_path
        self.make_cmd()
        return self.RunBlast(Verbose=Verbose)
        
        
#BR = BlastRunObj(query_path, '', 'nodb', dust="no", outfmt="7", num_threads="8", word_size="16", num_descriptions="10000")
#BR.Blast2Seqs(subj_path, out_path, params={})

