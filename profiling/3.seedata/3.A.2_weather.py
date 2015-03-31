
install.packages("weatherData")
library("weatherData")


WeatherData[{45, -90}, "WindSpeed"]
getWeatherForDate("KVIRGIN2", "2014-05-05")

head(seeds)


import pyowm

owm = pyowm.OWM('ed487e997b006bdb01e22853aef4883a')

# Will it be sunny tomorrow at this time in Milan (Italy) ?
forecast = owm.daily_forecast("Milan,it")
tomorrow = pyowm.timeutils.tomorrow()
forecast.will_be_sunny_at(tomorrow)  # Always True in Italy, right? ;-)

# Search for current weather in London (UK)
observation = owm.weather_at_place('London,uk')
w = observation.get_weather()
print(w)                      # <Weather - reference time=2013-12-18 09:20, 
# status=Clouds>

# Weather details
w.get_wind()                  # {'speed': 4.6, 'deg': 330}
w.get_humidity()              # 87
w.get_temperature('celsius')  # {'temp_max': 10.5, 'temp': 9.7, 'temp_min': 9.0}

# Search current weather observations in the surroundings of 
# lat=22.57W, lon=43.12S (Rio de Janeiro, BR)
observation_list = owm.weather_around_coords(-22.57, -43.12)

obs = owm.weather_at_coords(18.58333, -98.83333)

owm.weather_history_at_place(18.58333, -98.83333, '2013-09-13 16:46:40+00', '2013-09-13 19:16:40+00')




