all: test.c model.c scratch.c rigid_bodies.c common.c
	gcc test.c common.c -o test -lm -O3 -Wall -Wextra
	gcc model.c common.c -o model -lm -O3 -Wall -Wextra
	gcc scratch.c common.c -lm -O3 -Wall -Wextra
	gcc rigid_bodies.c common.c -o rigid -lm -O3 -Wall -Wextra

test: test.c common.c 
	gcc test.c common.c -o test -lm -O3 -Wall -Wextra
	
model: model.c common.c 
	gcc model.c common.c -o model -lm -O3 -Wall -Wextra
	
scratch: scratch.c common.c
	gcc scratch.c common.c -lm -O3 -Wall -Wextra
	
rigid: rigid_bodies.c common.c
	gcc rigid_bodies.c common.c -o rigid -lm -O3 -Wall -Wextra
