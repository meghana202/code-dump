using ArgParse
using CodecZlib
using FASTX
using BioSequences
using HTTP

cli = ArgParseSettings()
@add_arg_table cli begin
        "--file"
                help = "file"
        "--fasta"
                help = "fasta"
end
arg = parse_args(ARGS, cli)

#= READER
reader = nothing
if splitext(arg["fasta"])[2] == ".gz"
	reader = FASTA.Reader(GzipDecompressorStream(open(arg["fasta"])))
else
	reader = FASTA.Reader(open(arg["f"]))
end
=#

#MPRT, file
function protein_motif()
	file = read(open(arg["file"], "r"), String)
	file = split(file, "\n")
	#used length(file)-1 because the dataset comes with an empty line
	file = file[1:length(file)-1]
	website = "https://www.uniprot.org/uniprot/"
	for ID in file 
		r = HTTP.open("GET", website*ID*".fasta") do io
   			while !eof(io)
   				x = String(readavailable(io))
   				x = split(x, "\n")
   				x = join(x[2:length(x)-1])
   				locations = [ID]
   				for i = 1:length(x)-3
   					protein_seq = x[i:i+3]
   					if protein_seq[1] == 'N' && protein_seq[2] != 'P' && (protein_seq[3] == 'S'||protein_seq[3]=='T') &&protein_seq[4]!=4
   						push!(locations, string(i))
   					end
   				end 
   				if length(locations) > 1
   					println(locations[1])
   					for i = 2:length(locations) print("$(locations[i]) ") end 
   					println()
   				end
   			end
		end
	end
end
#protein_motif()

amino_acids = read(open("aa.txt", "r"), String)
amino_acids = split(amino_acids, ('\n',' '))
amino_acids = [x for x in amino_acids if x != ""]
aa = Dict()
for i = 2:2:length(amino_acids)
	if haskey(aa,amino_acids[i]) aa[(amino_acids[i])]+=1
	else aa[(amino_acids[i])] = 1 end
end
#println(aa)

#MRNA, uses aa.txt
function mRNA()
	count_ = 3
	file = read(open(arg["file"], "r"), String)
	file = join(split(file, '\n'))
	begin_ = findfirst(isequal('M'), file)
	end_ = findfirst(isequal("Stop"), file)
	if end_ != nothing
		protein = file[begin_:end_-1]
	elseif end_ == nothing
		protein = file[begin_:end]
	end
	for letter in protein
		count_ *= aa[string(letter)]
		mod_ = count_%1000000
		count_ = mod_
	end
	println(count_)
end
#mRNA()

function PERM(n)
	num_ = factorial(n)
	println(num_)
	answer = []
	if n == 1
		return [1]
	elseif n == 2
		return [1 2; 2 1;]
	elseif n > 2
		answer = zeros(factorial(n), n)
		for (i,j) in PERM(n-1)
			
			
		end
		return answer
	end
end
println(PERM(4))


