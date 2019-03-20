paraboloids.so: paraboloids.c
	gcc -Wall -fPIC -shared $< -o $@

.PHONY: clean
clean:
	rm -f paraboloids.so
