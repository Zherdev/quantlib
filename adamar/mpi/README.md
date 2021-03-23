Однокубитное преобразование Адамара. MPI.
============================================

```
$ git clone https://github.com/Zherdev/quantum.git
$ cd quantum/adamar/mpi
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
$ ./bin/adamar transform qubits_num target_qubit_num input_filename output_filename [-t]
```

Сгенерировать вектор состояния:

```
$ ./bin/adamar generate qubits_num output_filename [-t]
```

Сравнение векторов состояний:

```
$ ./bin/adamar cmp qubits_num a_filename b_filename
```

Флаг `-t` - для тестирования, вывод доп. информации в stdout.


[Отчет о запуске на машине IBM Polus](https://github.com/Zherdev/quantum/blob/master/adamar/mpi/report.pdf)
-------
