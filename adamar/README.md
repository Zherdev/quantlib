Однокубитное преобразование Адамара. OpenMP.
============================================

```
$ git clone https://github.com/Zherdev/quantum.git
$ cd quantum/adamar
```

Сборка:
```
$ make
```

Запуск:

```
$ ./bin/adamar N K thread_num [-t]
```

Флаг `-t` - для тестирования, вывод доп. информации в stdout.
`N` - число кубитов, `K` - номер выбранного кубита.

Пример:

```
$ ./bin/adamar 2 1 1 -t
```
```
Quantum Adamar transformation program started.
Random vector initialization started...
Random vector initialization is done.
Quant transform started...
Quant transform is done.
Time elapsed: 0.000004

Args:
        N: 2
        K: 1
        Thread num: 1
        Is test: 1

State vector v:
        Qubits num: 2
        Quantum state vector len: 4
        Float complex size: 8
        Quantum state vector size: 32
            0: (4.0000, 4.0000)
            1: (2.0000, 4.0000)
            2: (3.0000, 3.0000)
            3: (4.0000, 3.0000)

State vector w:
        Qubits num: 2
        Quantum state vector len: 4
        Float complex size: 8
        Quantum state vector size: 32
            0: (4.9497, 4.9497)
            1: (4.2426, 4.9497)
            2: (0.7071, 0.7071)
            3: (-1.4142, 0.7071)
```

[Отчет о запуске на машине Polus](https://github.com/Zherdev/quantum/blob/master/adamar/report.pdf)
-------
