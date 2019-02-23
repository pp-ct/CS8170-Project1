import os;

def write_to_file(file_name, data, mode):
	with open(file_name, mode) as file:  
		file.write(data+"\n") 

def write_stream(filename,data_stream):
	with open(filename, "w+") as handle:
	    handle.write(data_stream.read())
	    data_stream.close()

def create_dir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
