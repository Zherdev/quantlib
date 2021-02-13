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

Пример:

```
./bin/adamar 28 10 16 -t
```

Флаг `-t` - для тестирования, вывод доп. информации в stdout.
`N` - число кубитов, `K` - номер выбранного кубита.

[Отчет](https://github.com/Zherdev/quantum/blob/master/adamar/report.pdf)
-------