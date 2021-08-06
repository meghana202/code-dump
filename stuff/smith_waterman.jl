#match = 5, mismatch = -1, gap = -1
function sw(seq1, seq2)
	matrix_ = fill(0,length(seq2)+1,length(seq1)+1)
	trace_ = Matrix{Union{Char, Nothing}}(nothing,length(seq2)+1,length(seq1)+1)

	#first row and column
	for i = 1:length(seq1)
		trace_[1, i+1] = 'L'
	end
	for i = 1:length(seq2)
		trace_[i+1, 1] = 'U'
	end
	trace_[1,1] = '-' 
	
	#fill
	max_score = 0
	max_i = 0
	max_j = 0
	for i = 2:length(seq2)+1
		for j = 2:length(seq1)+1
			#x stores whether it was a match or mismatch/indel
			if seq1[j-1]==seq2[i-1] x = 5 else x = -1 end
			
			#filling matrix
			top = matrix_[i-1,j] - 1
			left = matrix_[i, j-1] - 1
			di = matrix_[i-1,j-1] + x
			
			if 0 > di && 0 > top && 0 > left 
				matrix_[i,j] = 0 
			elseif di > top && di > left
				if di > max_score 
					max_score = di
					max_i = i
					max_j = j
				end
				matrix_[i,j] = di
				trace_[i,j] = 'D'
			elseif top > left
				matrix_[i,j] = top
				trace_[i,j] = 'U'
			else 
				matrix_[i,j] = left
				trace_[i,j] = 'L'
			end			
		end
	end

	#trace
	seqs = []
	x = max_i
	y = max_j
	loc = findall(x->x == max_score, matrix_)
	length_ = length(loc)
	for i = 1:length_
		cell_ = loc[i]
		firstseq = ""
		secondseq = ""
		while true
			x = cell_[1]
			y = cell_[2]
			cell_ = [x, y]
			if matrix_[x,y] == 0 break end
			if trace_[x,y] == 'U'
				secondseq *= seq2[x-1]
				firstseq *= "-"
				cell_ = [x-1,y]
			elseif trace_[x,y] == 'L'
				firstseq *= seq1[y-1]
				secondseq *= "-"
				cell_ = [x, y-1]
			elseif trace_[x,y] == 'D'
				secondseq *= seq2[x-1]
				firstseq *= seq1[y-1]
				cell_ = [x-1, y-1]
			end
		end
		push!(seqs,(firstseq, secondseq))
	end
	
	#output/formatting
	longest = 0 
	for i = 1:length(seqs)
		if length(seqs[i][1]) > longest
			longest = length(seqs[i][1])
			global answer = seqs[i]
		end
	end
	for element in answer println(reverse(element)) end
	
end
sw(ARGS[1], ARGS[2])
