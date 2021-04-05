Однокубитное квантовое преобразование n-Адамар с зашумленными вентилями.
MPI + OpenMP.
========================================================================

```
$ git clone https://github.com/Zherdev/quantum.git
$ cd quantum/hadamard/mpi
```

Сборка:
```
$ make
```

Генерация данных:
```
$ make generate
```

Запуск для замеров времени работы алогоритма на сгенерированных данных:
```
$ make time
```

Запуск для замеров потери точности:
```
$ make loss
```

Выполнить преобразование c шумами:

```
$ ./bin/hadamard transform qubits_num threads_num accuracy input_filename output_filename [-t] [-l]
```

Сгенерировать вектор состояния:

```
$ ./bin/hadamard generate qubits_num threads_num output_filename [-t]
```


Флаг `-t` - для тестирования, вывод доп. информации в stdout.
Флаг `-l` - требуется ли подсчет потери точности.


[Отчет о запуске на машине Bluegene HPC](https://github.com/Zherdev/quantum/blob/master/hadamard_noise/report.pdf)
-------
