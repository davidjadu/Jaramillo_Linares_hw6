all: trayectoria_10.0_30.0_runge_kutta3d.pdf trayectoria_10.0_30.0_runge_kutta.pdf

trayectoria_10.0_30.0_runge_kutta3d.pdf: trayectoria_10.0_30.0_runge_kutta.pdf

trayectoria_10.0_30.0_runge_kutta.pdf: trayectoria_10.0_30.0_runge_kutta.dat plot_trayectoria.py
	python plot_trayectoria.py trayectoria_10.0_30.0_runge_kutta.dat

trayectoria_10.0_30.0_runge_kutta.dat: particle_in_field.x
	./particle_in_field.x 10 30

particle_in_field.x: particle_in_field.c
	cc particle_in_field.c -o particle_in_field.x -lm

clear:
	rm -f trayectoria*
	rm -f particle_in_field.x
