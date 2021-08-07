import itertools 
don = [50, 226, 345, 478]
acc = [120, 236, 378, 668]

for n in range(1,5):
	for dsites in itertools.combinations(don, n):
		for asites in itertools.combinations(acc,n):
			for i in range(n):
				
				print(dsites[i], asites[i], end=", ") 
			
			print()