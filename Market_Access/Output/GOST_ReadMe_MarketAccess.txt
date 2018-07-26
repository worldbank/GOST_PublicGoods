
        GOST: Market Access: Product Assumptions

        Last Updated: 26 Jul 2018
        Programmer: C. Fox
        Theory: K. Garrett, T. Norman

        This GOST Market Access product is based off of:
                - Mapbox's Matrix API for travel times;
                - OSRM's API for travel times

        Travel Time Calculation

        The Mapbox Matrix API provides estimated trip durations in seconds.
        The time it takes to travel from one point to another is determined by a
        number of factors, including:
        - The profile used (walking, cycling, or driving); (GOST: set to driving)
        - The speed stored in the maxspeed tag in OpenStreetMap
          (https://wiki.openstreetmap.org/wiki/Key:maxspeed)
        - Traffic derived from real-time telemetry data, provided by Mapbox

        Traffic data

        In addition to the contributions of OpenStreetMap, Mapbox SDKs collect
        anonymous data, or telemetry, about devices using their services to continuously
        update their routing network. Attributes such as speed, turn restrictions, and
        travel mode can be collected to improve OpenStreetMap.

        Advanced - Speed Assumptions

        See https://github.com/Project-OSRM/osrm-backend/blob/master/docs/profiles.md
        For a full explanation of profiles, and how speeds are calculated across segments

        Note on API request timings

        Requests using mapbox/driving, mapbox/walking, and mapbox/cycling profiles
        can specify up to 25 input coordinates per request. Requests using the
        mapbox/driving-traffic profiles can specify up to 10 input coordinates per request.

        Requests using mapbox/driving, mapbox/walking, and mapbox/cycling profiles
        have a maximum limit of 60 requests per minute. Requests using the
        mapbox/driving-traffic profiles have a maximum of 30 requests per minute.

        Algorithm flags

        Commands recognised for this script:
        -p   Path - string for input and output folder path
        -f   File name of .csv containing input data
        -m   Latitude column name.
        -n   Longitude column name
        -o   Origin Unique Identifier column name (e.g. District, Name, Object ID...).
             This is mainly helpful for joining the output back to the input data / a shapefile,
             and is non-essential in terms of the calculation. It can be text or a number.
        -q   Population / weighting column name
        -c   Server call type - "OSRM" for OSRM, "MB" for Mapbox, "MBT" for Mapbox traffic, or "Euclid" for Euclidean distances (as the crow flies)
        -l   Limit - use this to limit the coordinate input list (int). Optional.

        *** Optional - if sources =/= destinations. Note - Unique identifier and Pop column names must remain the same ***
        -W   Filename of destinations csv
        *** Optional - if encountering server errors / internet connectivity instability ***
        -R   Save - input latest save number to pick up matrix construciton process from there.
        -Z   Rescue number parameter - If you have already re-started the download process, denote how many times. First run = 0, restarted once = 1...
        Do NOT put column names or indeed any input inside quotation marks.
        The only exceptions is if the file paths have spaces in them.
        