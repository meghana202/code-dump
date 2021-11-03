import random


body = [" O \n", "/", "|", "\ \n", "/", " \ \n"]
words = ["zebra", "goat", "chicken"]

word = words[random.randint(0, len(words)-1)]
print(word)
answer = ' '.join(list(word)) + " "


for i in range(6):
	print(body[i], end = "")
	
		
x = 0
letters = ["_ "]*len(word)
print(''.join(letters))
wrong_guesses = []
while(x < 6):
	guess = input("Letter: ")
	if guess.isalpha == False:
		print("Invalid input!")
		continue
	if len(guess) != 1:
		print("Invalid input")
		continue
	if guess in wrong_guesses or guess in ''.join(letters):
		print("You already guessed that!")
		continue
	if guess in word:
		for i in range(len(word)):
			if word[i] == guess:
				letters[i] = guess + " "
		print(''.join(letters))
		for i in range(x):
			print(body[i], end = "")
		print()
	else:
		x += 1
		print(''.join(letters))
		for i in range(x):
			print(body[i], end = "")
		print()
		wrong_guesses.append(guess)
	print(f"Wrong Guesses: {' '.join(wrong_guesses)}")
	
		
	if answer == ''.join(letters):
		print("Nice job!")
		break

if x == 6: print("You died :(")