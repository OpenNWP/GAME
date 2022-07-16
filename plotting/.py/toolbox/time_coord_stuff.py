import datetime

def find_time_coord(year, month, day, hour, minute, second, microsecond):
	time_coord = datetime.datetime(year, month, day, hour, minute).replace(tzinfo = datetime.timezone.utc).timestamp()
	time_coord = time_coord + second + 1e-3*microsecond
	return time_coord

def return_date(time_coordinate):
	datetime_object = datetime.datetime.utcfromtimestamp(time_coordinate)
	year_found = datetime_object.year
	month_found = datetime_object.month
	day_found = datetime_object.day
	hrs_found = datetime_object.hour
	minute_found = datetime_object.minute
	second_found = datetime_object.second
	microsecond_found = datetime_object.microsecond
	return year_found, month_found, day_found, hrs_found, minute_found, second_found, microsecond_found
