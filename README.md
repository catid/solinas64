# Solinas64
## Finite field arithmetic modulo the Solinas prime p = 2^61-2^32+1 in C++

### Disclaimer: This is research software and not ready for production use.

The Solinas64 software takes advantage of the commonly available fast
64x64->128 multiplier to accelerate finite (base) field arithmetic.
So it runs a lot faster when built into a 64-bit executable.

I wrote this as a challenge to see if I could find a prime field that could run as fast as a Galois field.  This means I only needed to go far enough to implement a toy erasure code encoder.

As an erasure code, Solinas64 will expand the input data by a number of bytes (as shown during benchmarking).

## 64-bit Benchmarks

This is comparing the gf256 library to solinas64 library for a software erasure code application.
This is representative of the type of performance that can be achieved by using either library to implement a convolutional code.

Main result: Solinas64 is only 2-5 times slower than gf256.

Almost 98% of the execution time is in the inner loop performing the muladd operation on the data.

    Testing file size = 10 bytes
    N = 2 :  gf256_MBPS=434 Solinas64_MBPS=392 Solinas64_OutputBytes=16
    N = 4 :  gf256_MBPS=645 Solinas64_MBPS=579 Solinas64_OutputBytes=16
    N = 8 :  gf256_MBPS=941 Solinas64_MBPS=689 Solinas64_OutputBytes=16
    N = 16 :  gf256_MBPS=1066 Solinas64_MBPS=761 Solinas64_OutputBytes=16
    N = 32 :  gf256_MBPS=1212 Solinas64_MBPS=814 Solinas64_OutputBytes=16
    N = 64 :  gf256_MBPS=1233 Solinas64_MBPS=852 Solinas64_OutputBytes=16
    N = 128 :  gf256_MBPS=1263 Solinas64_MBPS=855 Solinas64_OutputBytes=16
    N = 256 :  gf256_MBPS=1286 Solinas64_MBPS=883 Solinas64_OutputBytes=16
    N = 512 :  gf256_MBPS=1282 Solinas64_MBPS=886 Solinas64_OutputBytes=16
    Testing file size = 100 bytes
    N = 2 :  gf256_MBPS=2222 Solinas64_MBPS=1904 Solinas64_OutputBytes=104.648
    N = 4 :  gf256_MBPS=3448 Solinas64_MBPS=2366 Solinas64_OutputBytes=105.216
    N = 8 :  gf256_MBPS=5161 Solinas64_MBPS=2523 Solinas64_OutputBytes=106.176
    N = 16 :  gf256_MBPS=6015 Solinas64_MBPS=2872 Solinas64_OutputBytes=107.64
    N = 32 :  gf256_MBPS=7785 Solinas64_MBPS=2777 Solinas64_OutputBytes=109.736
    N = 64 :  gf256_MBPS=8791 Solinas64_MBPS=2791 Solinas64_OutputBytes=111.456
    N = 128 :  gf256_MBPS=9509 Solinas64_MBPS=2774 Solinas64_OutputBytes=111.968
    N = 256 :  gf256_MBPS=9356 Solinas64_MBPS=2772 Solinas64_OutputBytes=112
    N = 512 :  gf256_MBPS=9458 Solinas64_MBPS=2634 Solinas64_OutputBytes=112
    Testing file size = 1000 bytes
    N = 2 :  gf256_MBPS=12738 Solinas64_MBPS=4338 Solinas64_OutputBytes=1000.62
    N = 4 :  gf256_MBPS=16806 Solinas64_MBPS=4228 Solinas64_OutputBytes=1001.26
    N = 8 :  gf256_MBPS=20100 Solinas64_MBPS=4284 Solinas64_OutputBytes=1002.41
    N = 16 :  gf256_MBPS=21857 Solinas64_MBPS=4286 Solinas64_OutputBytes=1003.69
    N = 32 :  gf256_MBPS=19900 Solinas64_MBPS=4288 Solinas64_OutputBytes=1005.93
    N = 64 :  gf256_MBPS=19518 Solinas64_MBPS=4273 Solinas64_OutputBytes=1007.39
    N = 128 :  gf256_MBPS=20340 Solinas64_MBPS=4190 Solinas64_OutputBytes=1007.98
    N = 256 :  gf256_MBPS=18030 Solinas64_MBPS=4102 Solinas64_OutputBytes=1008
    N = 512 :  gf256_MBPS=17460 Solinas64_MBPS=4066 Solinas64_OutputBytes=1008
    Testing file size = 10000 bytes
    N = 2 :  gf256_MBPS=25873 Solinas64_MBPS=4835 Solinas64_OutputBytes=10000.5
    N = 4 :  gf256_MBPS=23134 Solinas64_MBPS=4638 Solinas64_OutputBytes=10001.2
    N = 8 :  gf256_MBPS=23902 Solinas64_MBPS=4683 Solinas64_OutputBytes=10002.3
    N = 16 :  gf256_MBPS=24334 Solinas64_MBPS=4621 Solinas64_OutputBytes=10003.8
    N = 32 :  gf256_MBPS=23616 Solinas64_MBPS=4600 Solinas64_OutputBytes=10006.1
    N = 64 :  gf256_MBPS=23580 Solinas64_MBPS=4564 Solinas64_OutputBytes=10007.4
    N = 128 :  gf256_MBPS=23839 Solinas64_MBPS=4437 Solinas64_OutputBytes=10007.9
    N = 256 :  gf256_MBPS=23833 Solinas64_MBPS=4497 Solinas64_OutputBytes=10008
    N = 512 :  gf256_MBPS=23512 Solinas64_MBPS=4551 Solinas64_OutputBytes=10008
    Testing file size = 100000 bytes
    N = 2 :  gf256_MBPS=22864 Solinas64_MBPS=4630 Solinas64_OutputBytes=100001
    N = 4 :  gf256_MBPS=22637 Solinas64_MBPS=4523 Solinas64_OutputBytes=100001
    N = 8 :  gf256_MBPS=22243 Solinas64_MBPS=4424 Solinas64_OutputBytes=100002
    N = 16 :  gf256_MBPS=22471 Solinas64_MBPS=4470 Solinas64_OutputBytes=100004
    N = 32 :  gf256_MBPS=23604 Solinas64_MBPS=4805 Solinas64_OutputBytes=100006
    N = 64 :  gf256_MBPS=22448 Solinas64_MBPS=4789 Solinas64_OutputBytes=100007
    N = 128 :  gf256_MBPS=18957 Solinas64_MBPS=4728 Solinas64_OutputBytes=100008
    N = 256 :  gf256_MBPS=13059 Solinas64_MBPS=4775 Solinas64_OutputBytes=100008
    N = 512 :  gf256_MBPS=10815 Solinas64_MBPS=4723 Solinas64_OutputBytes=100008

Note that near the end it looks like the file sizes are exceeding the processor cache and it starts slowing down by 2x.

## 32-bit Benchmarks

When 64-bit operations are not available it is about 4x slower, which makes sense because it needs to use 4 multiplies instead of 1.

    Testing file size = 10 bytes
    N = 2 :  gf256_MBPS=183 Solinas64_MBPS=142 Solinas64_OutputBytes=16
    N = 4 :  gf256_MBPS=264 Solinas64_MBPS=192 Solinas64_OutputBytes=16
    N = 8 :  gf256_MBPS=281 Solinas64_MBPS=253 Solinas64_OutputBytes=16
    N = 16 :  gf256_MBPS=374 Solinas64_MBPS=296 Solinas64_OutputBytes=16
    N = 32 :  gf256_MBPS=438 Solinas64_MBPS=323 Solinas64_OutputBytes=16
    N = 64 :  gf256_MBPS=396 Solinas64_MBPS=335 Solinas64_OutputBytes=16
    N = 128 :  gf256_MBPS=440 Solinas64_MBPS=339 Solinas64_OutputBytes=16
    N = 256 :  gf256_MBPS=432 Solinas64_MBPS=336 Solinas64_OutputBytes=16
    N = 512 :  gf256_MBPS=436 Solinas64_MBPS=342 Solinas64_OutputBytes=16
    Testing file size = 100 bytes
    N = 2 :  gf256_MBPS=2222 Solinas64_MBPS=673 Solinas64_OutputBytes=104.648
    N = 4 :  gf256_MBPS=2721 Solinas64_MBPS=760 Solinas64_OutputBytes=105.216
    N = 8 :  gf256_MBPS=4232 Solinas64_MBPS=856 Solinas64_OutputBytes=106.176
    N = 16 :  gf256_MBPS=4610 Solinas64_MBPS=855 Solinas64_OutputBytes=107.64
    N = 32 :  gf256_MBPS=4892 Solinas64_MBPS=896 Solinas64_OutputBytes=109.736
    N = 64 :  gf256_MBPS=5470 Solinas64_MBPS=905 Solinas64_OutputBytes=111.456
    N = 128 :  gf256_MBPS=5643 Solinas64_MBPS=920 Solinas64_OutputBytes=111.968
    N = 256 :  gf256_MBPS=5693 Solinas64_MBPS=916 Solinas64_OutputBytes=112
    N = 512 :  gf256_MBPS=5798 Solinas64_MBPS=916 Solinas64_OutputBytes=112
    Testing file size = 1000 bytes
    N = 2 :  gf256_MBPS=10152 Solinas64_MBPS=1159 Solinas64_OutputBytes=1000.62
    N = 4 :  gf256_MBPS=12461 Solinas64_MBPS=1194 Solinas64_OutputBytes=1001.26
    N = 8 :  gf256_MBPS=14388 Solinas64_MBPS=1145 Solinas64_OutputBytes=1002.41
    N = 16 :  gf256_MBPS=15794 Solinas64_MBPS=1156 Solinas64_OutputBytes=1003.69
    N = 32 :  gf256_MBPS=15488 Solinas64_MBPS=1167 Solinas64_OutputBytes=1005.93
    N = 64 :  gf256_MBPS=15355 Solinas64_MBPS=1175 Solinas64_OutputBytes=1007.39
    N = 128 :  gf256_MBPS=15720 Solinas64_MBPS=1179 Solinas64_OutputBytes=1007.98
    N = 256 :  gf256_MBPS=14573 Solinas64_MBPS=1175 Solinas64_OutputBytes=1008
    N = 512 :  gf256_MBPS=14352 Solinas64_MBPS=1181 Solinas64_OutputBytes=1008
    Testing file size = 10000 bytes
    N = 2 :  gf256_MBPS=23446 Solinas64_MBPS=1275 Solinas64_OutputBytes=10000.5
    N = 4 :  gf256_MBPS=24154 Solinas64_MBPS=1246 Solinas64_OutputBytes=10001.2
    N = 8 :  gf256_MBPS=24539 Solinas64_MBPS=1227 Solinas64_OutputBytes=10002.3
    N = 16 :  gf256_MBPS=24416 Solinas64_MBPS=1231 Solinas64_OutputBytes=10003.8
    N = 32 :  gf256_MBPS=23400 Solinas64_MBPS=1230 Solinas64_OutputBytes=10006.1
    N = 64 :  gf256_MBPS=23615 Solinas64_MBPS=1233 Solinas64_OutputBytes=10007.4
    N = 128 :  gf256_MBPS=23829 Solinas64_MBPS=1272 Solinas64_OutputBytes=10007.9
    N = 256 :  gf256_MBPS=24241 Solinas64_MBPS=1303 Solinas64_OutputBytes=10008
    N = 512 :  gf256_MBPS=24110 Solinas64_MBPS=1306 Solinas64_OutputBytes=10008
    Testing file size = 100000 bytes
    N = 2 :  gf256_MBPS=23148 Solinas64_MBPS=1317 Solinas64_OutputBytes=100001
    N = 4 :  gf256_MBPS=23772 Solinas64_MBPS=1360 Solinas64_OutputBytes=100001
    N = 8 :  gf256_MBPS=23432 Solinas64_MBPS=1319 Solinas64_OutputBytes=100002
    N = 16 :  gf256_MBPS=23689 Solinas64_MBPS=1324 Solinas64_OutputBytes=100004
    N = 32 :  gf256_MBPS=23356 Solinas64_MBPS=1316 Solinas64_OutputBytes=100006
    N = 64 :  gf256_MBPS=21559 Solinas64_MBPS=1313 Solinas64_OutputBytes=100007
    N = 128 :  gf256_MBPS=18117 Solinas64_MBPS=1312 Solinas64_OutputBytes=100008
    N = 256 :  gf256_MBPS=12393 Solinas64_MBPS=1309 Solinas64_OutputBytes=100008
    N = 512 :  gf256_MBPS=10855 Solinas64_MBPS=1308 Solinas64_OutputBytes=100008

## API

Supported arithmetic operations: Add, Subtract, Multiply, Mul Inverse.  See solinas64.h.

There are also bulk memory operations useful for erasure codes: MultiplyRegion, MultiplyAddRegion.


## Future Work

TODO:
+ Write unit tests and validate all the arithmetic operations.
+ Implement an erasure code decoder to allow for end-to-end validation and real use.


#### Credits

Software by Christopher A. Taylor <mrcatid@gmail.com>.

Please reach out if you need support or would like to collaborate on a project.
