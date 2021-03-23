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
Time elapsed: 0.000005

Args:
        N: 2
        K: 1
        Thread num: 1
        Is test: 1

State vector v:
        Qubits num: 2
        Quantum state vector len: 4
        double complex size: 8
        Quantum state vector size: 32
            0: (0.4104, 0.4104)
            1: (0.2052, 0.4104)
            2: (0.3078, 0.3078)
            3: (0.4104, 0.3078)

State vector w:
        Qubits num: 2
        Quantum state vector len: 4
        double complex size: 8
        Quantum state vector size: 32
            0: (0.5078, 0.5078)
            1: (0.4353, 0.5078)
            2: (0.0725, 0.0725)
            3: (-0.1451, 0.0725)
```

[Отчет о запуске на машине IBM Polus](https://github.com/Zherdev/quantum/blob/master/adamar/omp/report.pdf)
-------
