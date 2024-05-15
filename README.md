# Autonomous-Robotics_Ex0 - GNSS Raw Mesurments:

## How to run:

* Go to --> filtered_csv\gnss-analysis-main\gnssutils 
* Enter filter.py (current directory)
* Insert your desired log file path into the constant "INPUT_LOG_FILE".
* Run filter.py

## Expected output:

(i) CSV file with the following additional columns: GPS time, SatPRN (ID), Sat.X, Sat.Y, Sat.Z, Pseudo-Range,
CN0, Doppler, Pos.X, Pos.Y, Pos,Z, Lat, Lon, Alt.

(ii) KML file with a computed path that passes between all the user's locations.
 

## Background

This assignment focuses on the basic principles of GNSS.
In particular, in this assignment an algorithm is implemented that receives parameters from a log file, 
such as satellite transmission time, satellite location, pseudo-range, etc. , 
and based on this data calculates the location of the user with the transmitter (could be a phone, or another technological device).

## Students

Amit Steinmetz - 207279373

Liron Cohen - 312324247

Bat-Ya Ashkenazi - 319088381

Maya Hadad - 209963784
