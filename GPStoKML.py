""" Author: AbuBakr Ghaznavi
    Description: Visualize the GPS data using KML
    Data: 11/29/2020
"""
# Import necessary files
import sys
from statistics import median
import numpy as np
import pandas as pd
import datetime
import re

# Define a regex to find multiples on a line
duplicate_finder = re.compile("GP.{3}")

# Class for giving generic data access for the command sentences
class GPCMD():
	# Get the longitude of the data
	def long(self):
		coeff = 1
		if self.horizontal == "W":
			coeff = -1
		return coeff * longFrac(self.longitude)
	# Get the latitude of the data
	def lat(self):
		coeff = 1
		if self.vertical == "S":
			coeff = -1	
		return coeff * latFrac(self.latitude)
	# Get the speed of the data
	def spd(self):
		return float(self.speed)
	# Only using the timestamp as a point of reference
	def ts(self):
		dt = datetime.datetime.strptime(self.timestamp, "%H%M%S.%f")
		return dt
	# Get the altitude of the data
	def alt(self):
		return float(self.antenna_height_geoid)
	# Get the angle of the data
	def angle(self):
		return float(self.true_course)
	# Get the number of the satellites
	def num_sat(self):
		return int(self.num_satellites)
# Get the GPGGA data
class GPGGA(GPCMD):
	def __init__(self, line):
		pieces = line.split(",")
		self.code = pieces[0].strip()
		self.timestamp = pieces[1].strip()
		self.latitude = pieces[2].strip()
		self.vertical = pieces[3].strip()
		self.longitude = pieces[4].strip()
		self.horizontal = pieces[5].strip()
		self.quality_indicator = pieces[6].strip()
		self.num_satellites = pieces[7].strip()
		self.hori_dilution = pieces[8].strip()
		self.antenna_height_geoid = pieces[9].strip()
		self.antenna_height_meters = pieces[10].strip()
		self.geoidal_seperation = pieces[11].strip()
		self.unit_seperation = pieces[12].strip()
		self.update_age = pieces[13].strip()
		self.diff_refr_id = pieces[14].strip()

	def __repr__(self):
		gps_map = {"Code": self.code,
			   "Timestamp": self.timestamp,
 			   "Latitude": self.latitude,
			   "Vertical": self.vertical,
			   "Longitude": self.longitude,
			   "Horizontal": self.horizontal,
			   "Number of satellites": self.num_satellites,
			   "Horizontal dillution": self.hori_dilution,
			   "Antenna height geoid": self.antenna_height_geoid,
			   "Antenna height meters": self.antenna_height_meters,
			   "Geoidal seperation unit":self.geoidal_seperation,
			   "Meters seperation": self.unit_seperation,
			   "Update age":self.update_age,
			   "Difference reference ID":self.diff_refr_id}
		string_format = ""
		for key in sorted(gps_map.keys()):
			string_format += key + ":" +  gps_map[key]
			string_format += "\n"	
		return string_format

# Get the Recommended minimum requirements
class GPRMC(GPCMD):
	def __init__(self, line):
		pieces = line.split(",")
		self.code = pieces[0].strip()
		self.timestamp = pieces[1].strip()
		self.validity = pieces[2].strip()
		self.latitude = pieces[3].strip()
		self.vertical = pieces[4].strip()
		self.longitude = pieces[5].strip()
		self.horizontal = pieces[6].strip()
		# Speed in knots
		self.speed = pieces[7].strip()
		self.true_course = pieces[8].strip()
		self.datestamp = pieces[9].strip()
		self.variation = None
		self.variation_horizontal = None
		if len(pieces) >  10:
			self.variation = pieces[10].strip()
		if len(pieces) > 11:
			self.variation_horizontal = pieces[11].strip()
	def __repr__(self):
		gps_map = {"Timestamp": self.timestamp,
			   "Validity": self.validity,	
			   "Latitude": self.latitude,	
			   "Vertical": self.vertical,	
			   "Longitude": self.longitude,	
			   "Horizontal": self.horizontal,	
			   "Speed": self.speed,	
			   "True course": self.true_course,	
			   "Datestamp": self.true_course,
			   "Variation": self.variation,
			   "Variation horizontal": self.variation_horizontal,
			   "Code": self.code}
		          
		string_format = ""
		for key in sorted(gps_map.keys()):
			val = gps_map[key]
			if val == None:
				str_val = "None"
			else:
				str_val = str(val)
			string_format += key + ":" + str_val
			string_format += "\n"	
		return string_format


# Define a method for reading the GPS recorder
class GPSReader():
	sentence_codes = ["$GPRMC","$GPGGA"]
	def __init__(self, filename):
		# Split the header into fields to read
		try:
			self.filehandle = open(filename, "r")
			self.version = self.filehandle.readline().replace("Vers","").strip()
			self.use_serial_feedback = self.filehandle.readline().split("=")[1].strip().lower() == "true"
			self.development_mode = self.filehandle.readline().split("=")[1].strip().lower() == "true"
			self.rmc_only = self.filehandle.readline().split("=")[1].strip().lower() == "true"
			# Read the blank line to align with the rest of the file
			self.filehandle.readline()
		except:
			print("Error: file is not valid GPS format")
			sys.exit(1)
	
	# Get the code from the line
	@staticmethod
	def getCode(line):
		comma_index = line.find(",")
		return line[0:comma_index]

	# Get the command sentence
	def getGPCmd(self):
		code = ""
		while code not in GPSReader.sentence_codes:
			
			line = self.filehandle.readline()
			# Skip over lines with double sentences
			if len(re.findall(duplicate_finder, line)) > 1:
				continue
			# If we reached EOF
			if len(line) == 0:
				return None
			code = GPSReader.getCode(line).strip()
				
		# Throw error if code is not GPS code
		if not(code in GPSReader.sentence_codes):
			raise Exception("Invalid GP code")
		if code == "$GPGGA":
			return GPGGA(line)
		elif code == "$GPRMC":
			return GPRMC(line)
	
	# Define an iterator to reading the lines as commands
	def iterCmds(self):	
		while True:
			GPCmd = self.getGPCmd()
			if GPCmd == None:
				self.filehandle.close()
				break
			yield GPCmd	

# Function to convert to fractional degrees
def longFrac(long_str):
	# Get first three digits for degree part
	deg_part = long_str[0:3]
	minute_part = long_str[3:]
	return float(deg_part) + float(minute_part)/60
	
	
# Convert the latitude to fractional	
def latFrac(long_str):
	# Get first three digits for degree part
	deg_part = long_str[0:2]
	minute_part = long_str[2:]
	return float(deg_part) + float(minute_part)/60
	
# Define a basic write to write the GPS data
class KMLWriter():
	HEADER = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
<Style id="yellowPoly">
  <LineStyle>
    <color>Afffff00</color>
    <width>6</width>
  </LineStyle>
  <PolyStyle>
    <color>7fffff00</color>
  </PolyStyle>
</Style>
<Placemark><styleUrl>#yellowPoly</styleUrl>
<LineString>
<Description>Speed in Knots, instead of altitude.</Description>
  <extrude>1</extrude>
  <tesselate>1</tesselate>
  <altitudeMode>absolute</altitudeMode>
  <coordinates>"""
	TRAILER = """
  </coordinates>
  </LineString>
  </Placemark>
 </Document>
</kml>
	"""
	def __init__(self,filename):
		self.filehandle = open(filename, "w")

	def writeFromDF(self, df):
		self.filehandle.write(KMLWriter.HEADER + "\n")
		for index,row in df.iterrows():
			write_vals = []
			for column_name in df.columns:
				write_vals.append(str(row[column_name]))
			write_str = ",".join(write_vals)	
			self.filehandle.write(write_str + "\n")
		self.filehandle.write(KMLWriter.TRAILER + "\n")
# Make a basic visualizer	


if len(sys.argv) > 2:
	TEST_FILENAME = sys.argv[1]
	OUT_FILENAME = sys.argv[2]
else:
	print("Usage GPStoKML.py <gps_file> <outfile>.kml")
	sys.exit(1)

SPEED_TOLERANCE = 0.1
ANGLE_SPEED_TOLERANCE = 0.01
# Have max speed tolerance to exclude lane changes
MAX_SPEED_TOLERANCE = 10

reader = GPSReader(TEST_FILENAME)

def angle_dist(a, b):
	d = a - b
	d = ( d + 180 ) % 360 - 180
	return d

def calcAngleDiffs(df):
	# Keep an initial zero for the sake of indexing
	look_ahead = 1
	diffs = [0] * look_ahead
	num_rows = len(df.index)
	for i in range(num_rows - look_ahead):
		row = df.iloc[i]
		next_row = df.iloc[i + look_ahead]
		first_angle = row["angle"]
		second_angle = next_row["angle"]
		# Calculate the angle delta
		angle_difference = angle_dist(second_angle, first_angle)
		# Apply the speed as a factor (average between them)
		# Get the average speed
		avg_speed = (row["speed"] + next_row["speed"]) / 2.0
		# Exclude angles when there is no speed factor
		if avg_speed < ANGLE_SPEED_TOLERANCE:
			angle_difference = 0
		if avg_speed > MAX_SPEED_TOLERANCE:
			angle_difference = 0		
		diffs.append(angle_difference)
	return np.array(diffs)

def find_stop_range(df):
	# Go to the end of the data frame and check the last time the speed is under the threshold
	last_index = len(df.index) - 1
	cur_idx = last_index
	for index, row in df[::-1].iterrows():
		current_speed = row["speed"]
		if current_speed > SPEED_TOLERANCE:
			return (last_index, cur_idx)
		cur_idx -= 1
	return (start_index, cur_idx + 1)

def find_start_range(df):
	start_index = 0
	cur_idx = start_index
	for index, row in df[::-1].iterrows():
		current_speed = row["speed"]
		if current_speed < SPEED_TOLERANCE:
			return (start_index, cur_idx)
		cur_idx += 1	
	# Leave one value in the range to represent
	return (start_index, cur_idx - 1)
	
def agglo_turns(df):
	# Use on the cross index of the dataframe and a sliding window approach for smoothing
	# Adjustable parameter for this function
	WIN_TOLERANCE = 3
	in_turn = False 
	window_size = 8
	middle_indices_arr = []
	num_rows = len(df.index)
	medians = []
	for idx in range(num_rows - window_size):
		middle_idx = int((idx + window_size) // 2)
		window = []
		for window_cell in range(window_size):
			cell_cross = df.iloc[idx + window_cell]["cross"]
			window.append(cell_cross)
		window_mean = np.mean(np.array(window))
		# Mean high enough to enter turn
		high_mean = window_mean > WIN_TOLERANCE
		# If mean was comprised of 1 high value surrouned by negatives then it will not count
		count_pos = len([num for num in window if num > 0])
		count_prop = count_pos / window_size
		if count_prop <= 0.7:
			high_mean = False	
		if not(in_turn) and high_mean:
			middle_indices_arr.append(middle_idx)
			in_turn = True
		if not(in_turn) and not(high_mean):
			# Do nothing if we aren't in high mean and not gonna be
			continue
		if in_turn and high_mean:
			# Append to our middle_indices
			middle_indices_arr.append(middle_idx)
		# Tough state where we handle transition to not in_turn
		# Save the median of our indices
		if in_turn and not(high_mean):
			in_turn = False
			medians.append(int(median(middle_indices_arr)))
			middle_indices_arr = []
	return medians		
		

def r_agglo_turns(df):
	# Use on the cross index of the dataframe and a sliding window approach for smoothing
	# Adjustable parameter for this function
	WIN_TOLERANCE = -3
	in_turn = False 
	window_size = 8
	middle_indices_arr = []
	num_rows = len(df.index)
	medians = []
	for idx in range(num_rows - window_size):
		middle_idx = int((idx + window_size) // 2)
		window = []
		for window_cell in range(window_size):
			cell_cross = df.iloc[idx + window_cell]["cross"]
			window.append(cell_cross)
		window_mean = np.mean(np.array(window))
		# Mean high enough to enter turn
		high_mean = window_mean < WIN_TOLERANCE
		# If mean was comprised of 1 high value surrouned by negatives then it will not count
		count_neg = len([num for num in window if num < 0])
		count_prop = count_neg / window_size
		if count_prop <= 0.7:
			high_mean = False	
		if not(in_turn) and high_mean:
			middle_indices_arr.append(middle_idx)
			in_turn = True
		if not(in_turn) and not(high_mean):
			# Do nothing if we aren't in high mean and not gonna be
			continue
		if in_turn and high_mean:
			# Append to our middle_indices
			middle_indices_arr.append(middle_idx)
		# Tough state where we handle transition to not in_turn
		# Save the median of our indices
		if in_turn and not(high_mean):
			in_turn = False
			medians.append(int(median(middle_indices_arr)))
			middle_indices_arr = []
	return medians		

data_mat = []
# Get all the latitudes and longitudes
cmds = [item for item in reader.iterCmds()]
for i in range(len(cmds) - 1):
	if i % 2 == 0:
		continue
	# Means that we have a pair of the same thing
	first_code = cmds[i].code
	second_code = cmds[i + 1].code
	if first_code == second_code:
		continue
	

	first_code = cmds[i].code
	second_code = cmds[i + 1].code

	

	if first_code == "$GPGGA":
		gpgga = cmds[i]
		gprmc = cmds[i + 1]
	elif first_code == "$GPRMC":
		gpgga = cmds[i + 1]
		gprmc = cmds[i]

	# Check for low fixquality or invalid 
	if gprmc.validity == "V" or gpgga.quality_indicator == "0":
		continue


	# Flatten out the relevant data fields to be used for pandas array
	# Most relevant fields will come from the GPRMC data
	data_row = [gprmc.ts(), gpgga.long(), gpgga.lat(), gpgga.alt(), gprmc.spd(), gprmc.angle(), gpgga.num_sat()]
	data_mat.append(data_row)

# Make the data frame from the data
column_list = ["time","longitude","latitude","altitude","speed","angle","num_sat"]
gps_data = pd.DataFrame(data_mat, columns=column_list)	
angle_dif = calcAngleDiffs(gps_data)
# Calculate angle differennces and keep the sign
gps_data["cross"] = 8 * np.square(angle_dif) * np.sign(angle_dif)

stop_range = find_stop_range(gps_data)
start_range = find_start_range(gps_data)
print(stop_range)
print(start_range)
# Drop non moving rows on both ends of the track
gps_data.drop(gps_data.index[start_range[0]:start_range[1]], inplace = True)
gps_data.drop(gps_data.index[stop_range[1]:stop_range[0]], inplace = True)

rs = agglo_turns(gps_data)
ls = r_agglo_turns(gps_data)

# Add zero column
gps_data["zero"] = np.zeros(len(gps_data.index))

print("Left:", rs)
print("Right:", ls)

writer = KMLWriter(OUT_FILENAME)
writer.writeFromDF(gps_data[["longitude","latitude","cross"]])
