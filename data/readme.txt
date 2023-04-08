Accelerometer, Gyroscope, Force Sensing Resistor Data (ADC) from Paddling Session on March 24, 2023
===================================================================================================
William Zhang
---------------------------------------------------------------------------------------------------
IMU = MPU6050
Sampling @ 100Hz
12bit ADC on STM32F401CUU6 Channel 1
---------------------------------------------------------------------------------------------------
There are 8 colummns:
1. Time			Time in ms since microcontroller started
2. Accelerometer X	Down the shaft of the paddle axis
3. Accelerometer Y	Left/Right of Paddle axis
4. Accelerometer Z	Forward/Back of Paddle axis
5. ADC			10kOhms voltage divider
6. Gyroscope X		Down the shaft of the paddle axis
7. Gyroscope Y		Left/Right of Paddle axis
8. Gyroscope Z		Forward/Back of Paddle axis
---------------------------------------------------------------------------------------------------
The raw data is contained in data.csv.
Workout was 8x(30s@5, 2min@7, 30s@5, 15s@9, 30s@5)
---------------------------------------------------------------------------------------------------
Data and Notes in Chronological Order:
1. mar24_2023_warmup.csv 
	- no wind
	- calmn water
	- light paddle
2. mar24_2023_set_1.csv
	- 30s@5, 2min@7, 30s@5, 15s@9, 30s@5
	- might had wrong hand placement
	- not the best piece in the world

3. mar24_2023_set_2.csv
	- 30s@5, 2min@7, 30s@5, 15s@9, 30s@5
	- ok peice
	- feel like IMU unit got shifted

4. mar24_2023_set_3.csv
	- 30s@5, 2min@7, 30s@5, 15s@9, 30s@5
	- not a bad set
	- broke the strap while trying to tighten it after this piece
	- slight tailwind
	- not orthgonal to blade

5. mar24_2023_set_4.csv
	- 30s@5, 2min@7, 30s@5, 15s@9, 30s@5
	- slight tailwind

6. mar24_2023_set_5.csv
	- 30s@5, 2min@7, 30s@5, 15s@9, 30s@5
	- moderate tailwind
	- recorded 

7. mar24_2023_set_6.csv
	- 30s@5, 2min@7, 30s@5, 15s@9, 30s@5
	- moderate headwind
	- recorded 

8. mar24_2023_set_7.csv
	- 30s@5, 2min@7, 30s@5, 15s@9, 30s@5
	- moderate headwind
	- good piece
	- paddled during the rest between set 7 adn set 8

9. mar24_2023_set_8.csv
	- 30s@5, 2min@7, 30s@5, 15s@9, 30s@5
	- moderate headwind
	- good piece
	- paddled during the rest between set 7 adn set 8

10. mar24_2023_normal.csv
	- normal paddling at PES 6
	- moderate headwind
	- choppy waves

11. mar24_2023_back.csv
	- flipped paddle backwards
	- FSR now on the back palm area

12. mar24_2023_left.csv
	- paddled on left side
	- better tech, less power?

13. mar24_2023_bad.csv
	- simulated bad paddling
	- blade not fully buried before pulling
	- werid motions
---------------------------------------------------------------------------------------------------