def read_data(file_name):
	file = []
	with open(file_name, "r") as filestream:
		for line in filestream:
			file.append([x.strip() for x in line.split(',')])
	
	#remove "time"
	file[0].pop(0) 
	file[1].pop(0)
	
	results={}
	for elem in range(len(file[0])):
		node, type = file[0][elem].split('.')
		if not node in results:
			results[node] = 0.0 + 0.0j
		if type=="real":
			results[node] += float(file[1][elem])
		elif type=="imag":
			results[node] += complex(0,float(file[1][elem]))
	
	return results