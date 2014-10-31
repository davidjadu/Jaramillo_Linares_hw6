all: trayectoria_10.0_30.0_runge_kutta3d.pdf trayectoria_10.0_30.0_runge_kutta.pdf

trayectoria_10.0_30.0_runge_kutta3d.pdf: trayectoria_10.0_30.0_runge_kutta.pdf

trayectoria_10.0_30.0_runge_kutta.pdf: trayectoria_10.0_30.0_runge_kutta.dat plot_trayectoria.py
	python plot_trayectoria.py trayectoria_10.0_30.0_runge_kutta.dat

trayectoria_10.0_30.0_runge_kutta.dat: a.out
	./a.out 10 30

a.out: particle_in_field.c
	cc particle_in_field.c

clear:
	rm -f trayectoria*
	rm -f a.out
