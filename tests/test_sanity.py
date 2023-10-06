import sanclone


def test_version():
    assert sanclone.__version__


def test_echo_tool():
    from sanclone import State
    from sanclone.tools import EchoTool

    tool = EchoTool(shared_state=State())
    assert tool.run("Hello") == "Hello"

# def test_state_tool():
#     from sanclone import State
#     from sanclone.State import download_genbank_file
#     accession_id_vector = "NC_005213"
#     output_filename_vector = "NC_005213.gbk"
#     accession_id_linear_insert = "NC_000932"
#     output_filename_linear_insert = "NC_000932.gbk"
#     download_genbank_file(accession_id_vector, output_filename_vector)
#     download_genbank_file(accession_id_linear_insert, output_filename_linear_insert)

#     for gb_record in SeqIO.parse(open(output_filename_linear_insert,"r"), "genbank") :
#         # now do something with the record
#         print ("Name %s, %i features" % (gb_record.name, len(gb_record.features)))

#     vector_seq = list(SeqIO.parse(open(output_filename_vector,"r"), "genbank"))
#     insert_seq = list(SeqIO.parse(open(output_filename_linear_insert,"r"), "genbank"))

#     seq_anno = State(vector_seq[0])
#     seq_anno.store_linear_insert(insert_seq[0])

#     retrieved_vector = seq_anno.retrieve_vector()
#     retrieved_insert = seq_anno.retrieve_linear_insert()
#     print(retrieved_vector)
#     print(retrieved_insert)