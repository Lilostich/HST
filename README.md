# MyVector

Учебный проект по предмету "***Гибридные суперкомпьютерные технологии***"

## Бизнес логика

**Расчет собственного вектора матрицы**

## Алгоритм

Расчет собственного вектора происходит с помощью [степенного метода](https://ru.wikipedia.org/wiki/%D0%A1%D1%82%D0%B5%D0%BF%D0%B5%D0%BD%D0%BD%D0%BE%D0%B9_%D0%BC%D0%B5%D1%82%D0%BE%D0%B4). Результатом является собственный вектор соответствующий максимальному собственному значению.

Основная формула: $ r_{k+1} = \frac{Ar_k}{||Ar_k||} $

## Вариант исполнения

* Последовательная реализация

## Особенности

1. Считывание данных происходит из файла (либо данные передаются по протоколу TCP)

2. Данные генерируются утилитой, принимающей в качестве параметров размер данных для обработки в мегабайтах и имя файла (TCP хост-порт) куда будут выгружены данные

3. Программа выполняет бизнес-логику и записывает результат в выходной файл (отправляет данные на порт возврата результатов программы-генератора по TCP, сохранение файла с результатами осуществляет программа-генератор)

4. В конце файла с результатами сохраняется информация о времени выполнения вычислений и размере обработанных данных
