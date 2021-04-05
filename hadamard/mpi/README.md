Однокубитное преобразование Адамара. MPI.
============================================

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

Тестирование:

```
$ make test
```

Запуск для замеров времени работы алогоритма на сгенерированных данных:
```
$ make generate
$ make run
```

Выполнить преобразование:

```
$ ./bin/hadamard transform qubits_num target_qubit_num input_filename output_filename [-t]
```

Сгенерировать вектор состояния:

```
$ ./bin/hadamard generate qubits_num output_filename [-t]
```

Сравнение векторов состояний:

```
$ ./bin/hadamard cmp qubits_num a_filename b_filename
```

Флаг `-t` - для тестирования, вывод доп. информации в stdout.


[Отчет о запуске на машине IBM Polus](https://github.com/Zherdev/quantum/blob/master/hadamard/mpi/report.pdf)
-------
