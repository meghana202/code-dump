using ArgParse

cli = ArgParseSettings()
@add_arg_table cli begin
        "--file"
                help = "file"
        "--fasta"
                help = "fasta"
end
arg = parse_args(ARGS, cli)


function fasta_decoder(file)
	file = readlines(open(file, "r"))
	id = []
	seqs = []
	i = 1
	while i <= length(file)
		seq = ""
		if file[i][1] == '>' 
			push!(id, file[i])
		end
		i += 1
		while file[i][1] != '>'
			seq = seq*file[i]
			i += 1
			if i > length(file) break end
		end
		push!(seqs, seq)
	end
	return(id, seqs)
end
#println(fasta_decoder("test.fa"))

#LONG, uses fasta

function superstring()
	dna = fasta_decoder(arg["fasta"])[2]
	x = length(dna[1])÷2+1
	answers = []
	for element in dna println(length(element)) end
	for i = x:length(dna[1])-1
		graph = Dict{String, String}()
		for element in dna
			graph[element[begin:i]] = element[length(element)-i+1:end]
		end 
		answer = []
		overlap = 2i - length(dna[1])
		for (key,value) in graph 
			for (key1, value1) in graph
				if length(answer) == length(dna) break end
				if value == key1 push!(answer, key, key1) end 
			end 
		end
		println()
		if length(answer) > 1 
			push!(answer, graph[answer[end]])
			push!(answers, answer)
			for j = 1:length(answer)
				if j == 1 print(answer[j]) 
				elseif j!=1
					print(answer[j][i÷2+2:end]) 
				end
			end
		end
	end
	println()
end
superstring()




























#REVP


function rc(dna)
	rc_ = ""
	for i = 0:length(dna)-1
		if dna[length(dna)-i] == 'A' rc_ *= 'T'
		elseif dna[length(dna)-i] == 'T' rc_ *= 'A'
		elseif dna[length(dna)-i] == 'G' rc_ *= 'C'
		elseif dna[length(dna)-i] == 'C' rc_ *= 'G' end
	end
	return rc_
end 
#rc("GAT")

function REVP()
	file = read(open(arg["file"], "r"), String)
	file = join(split(file,'\n')[2:end])
	for i = 1:length(file)-11
		for j = 4:2:12
			if file[i:i+j-1] == rc(file[i:i+j-1]) println("$i $j") end
		end
	end
	last_ = file[length(file)-10:end]
	for i = 1:11
		for j = 4:2:12-i
			if last_[i:i+j-1] == rc(last_[i:i+j-1]) println("$(length(file)-(11-i)) $j") end
		end
	end
		
end
#REVP()



function fasta_decoder(file)
	file = readlines(open(file, "r"))
	id = []
	seqs = []
	i = 1
	while i <= length(file)
		seq = ""
		if file[i][1] == '>' 
			push!(id, file[i])
		end
		i += 1
		while file[i][1] != '>'
			seq = seq*file[i]
			i += 1
			if i > length(file) break end
		end
		push!(seqs, seq)
	end
	return(id, seqs)
end
#println(fasta_decoder("test.fa"))




#SPLC
amino = Dict{Any, Any}("GAC" => "D", 
"UAU" => "Y", "UUA" => "L", "CUU" => "L", "CCU" => "P", 
"GGC" => "G", "GCU" => "A", "CGG" => "R", "CAG" => "Q", 
"ACG" => "T", "GAU" => "D", "UGU" => "C", "AGA" => "R", 
"GUU" => "V", "CUA" => "L", "UUG" => "L", "CUC" => "L", 
"CGA" => "R", "CAA" => "Q", "UUU" => "F", "UUC" => "F", 
"AUU" => "I", "CUG" => "L", "UAA" => "Stop", "AAC" => "N", 
"GUC" => "V", "UCC" => "S", "CCA" => "P", "CAC" => "H", 
"GCA" => "A", "CCG" => "P", "GAA" => "E", "GUA" => "V", 
"ACA" => "T", "GUG" => "V", "UCU" => "S", "AGG" => "R", 
"AUC" => "I", "AUG" => "M", "AUA" => "I", "UCG" => "S", 
"AGC" => "S", "ACU" => "T", "GCC" => "A", "AAG" => "K", 
"GGU" => "G", "AAU" => "N", "AAA" => "K", "UAG" => "Stop", 
"UGC" => "C", "ACC" => "T", "UGA" => "Stop", "AGU" => "S", 
"CCC" => "P", "UAC" => "Y", "GAG" => "E", "CGU" => "R", 
"GCG" => "A", "CAU" => "H", "CGC" => "R", "UGG" => "W",
"GGA" => "G", "UCA" => "S", "GGG" => "G")
function splicing()
	file = fasta_decoder(arg["file"])
	dna = file[2][1]
	introns = file[2][2:end]
	rna = ""
	for intron in introns
		dna = replace(dna, intron => "")
	end
	for letter in dna
		if letter == 'T' rna = rna*'U'
		else rna = rna*letter
		end
	end
	for i = 1:3:length(rna)-2
		if amino[rna[i:i+2]] == "Stop" break end
		print(amino[rna[i:i+2]])
	end
	println()
end
#splicing()
		