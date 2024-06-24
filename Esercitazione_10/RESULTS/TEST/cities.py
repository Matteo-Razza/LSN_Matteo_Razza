from geopy.geocoders import Nominatim
import time

def get_city_coordinates(city_names):
    geolocator = Nominatim(user_agent="city_coordinates")
    coordinates = {}
    
    for city in city_names:
        location = geolocator.geocode(city)
        if location:
            coordinates[city] = (location.latitude, location.longitude)
        else:
            coordinates[city] = None
        time.sleep(1)  # Sleep to avoid overwhelming the geocoding service
    
    return coordinates

city_names = ["Salerno, Italy", "Medio Campidano, Italy", "Sassari, Italy", "Savona, Italy", "Siena, Italy", "Syracuse, Italy", "Sondrio, Italy", "Taranto, Italy", "Teramo, Italy", "Terni, Italy", "Turin, Italy", "Ogliastra, Italy", "Trapani, Italy", "Trento, Italy", "Treviso, Italy", "Trieste, Italy", "Udine, Italy", "Varese, Italy", "Venice, Italy", "Verbania, Italy", "Vercelli, Italy", "Verona, Italy", "Vibo Valentia, Italy", "Vicenza, Italy", "Viterbo, Italy"]
coordinates = get_city_coordinates(city_names)

for city, coord in coordinates.items():
    print(f"Coordinates of {city}: {coord}")

#["Agrigento, Italy", "Alessandria, Italy", "Ancona, Italy", "Aosta, Italy", "Arezzo, Italy", "Ascoli Piceno, Italy", "Asti, Italy", "Avellino, Italy", "Bari, Italy", "Barletta, Italy", "Belluno, Italy", "Benevento, Italy", "Bergamo, Italy", "Biella, Italy", "Bologna, Italy", "Bolzano, Italy", "Brescia, Italy", "Brindisi, Italy", "Cagliari, Italy", "Caltanissetta, Italy", "Campobasso, Italy", "Iglesias, Italy", "Caserta, Italy", "Catania, Italy", "Catanzaro, Italy", "Chieti, Italy", "Como, Italy", "Cosenza, Italy", "Cremona, Italy", "Crotone, Italy", "Cuneo, Italy", "Enna, Italy", "Fermo, Italy", "Ferrara, Italy", "Florence, Italy", "Foggia, Italy", "Cesena, Italy", "Frosinone, Italy", "Genoa, Italy", "Gorizia, Italy", "Grosseto, Italy", "Imperia, Italy", "Isernia, Italy", "La Spezia, Italy", "L'Aquila, Italy", "Latina, Italy", "Lecce, Italy", "Lecco, Italy", "Livorno, Italy", "Lodi, Italy", "Lucca, Italy", "Macerata, Italy", "Mantua, Italy", "Carrara, Italy", "Matera, Italy", "Messina, Italy", "Milan, Italy", "Modena, Italy", "Monza, Italy", "Naples, Italy", "Novara, Italy", "Nuoro, Italy", "Olbia, Italy", "Oristano, Italy", "Padua, Italy", "Palermo, Italy", "Parma, Italy", "Pavia, Italy", "Perugia, Italy", "Pesaro, Italy", "Pescara, Italy", "Piacenza, Italy", "Pisa, Italy", "Pistoia, Italy", "Pordenone, Italy", "Potenza, Italy", "Prato, Italy", "Ragusa, Italy", "Ravenna, Italy", "Reggio Calabria, Italy", "Reggio Emilia, Italy", "Rieti, Italy", "Rimini, Italy", "Rome, Italy", "Rovigo, Italy"
