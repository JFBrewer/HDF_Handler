# Conda or pip installs
import pandas as pd
from pyhdf import SD
import numpy as np
from h5attr import H5Attr

# my packages 
from MDJ2KtoGreg import MDJ22K_Converter

class NDACC_Reader:
    def __init__(self):
        """
        Class constructor.
        """
        print("Calling __init__() constructor ...")
        

    def extract_fields_from_hdf4(self, file_name):
        """
        This method takes an HDF4 file and extracts all information from it
        It returns a dictionary with each key-variable pair
        """

        # Open the HDF4 file
        hdf_file = SD.SD(file_name)

        # Get all datasets in the file
        datasets = hdf_file.datasets()

        # Create an empty dictionary to store the data
        data_dict = {}

        # Loop through each dataset
        for dataset_name in datasets:
            # Select the current dataset
            dataset = hdf_file.select(dataset_name)
            # Get the data from the dataset
            dataset_data = np.array(dataset)
            # Add the data to the dictionary with the dataset name as the key
            data_dict[dataset_name] = dataset_data

        # Close the HDF4 file
        hdf_file.end()
    
        # Return the dictionary of data
        return data_dict

    def extract_fields_from_hdf5(self, file_name):
        """
        This method takes an HDF5 file and extracts all information from it
        It returns a dictionary with each key-variable pair
        It is analagous to the HDF4 routine and attempts to return the same kind of object
        """
        # Load in the HDF5 file - must specify not lazy loading
        hdf5_file = H5Attr(file_name, lazy = False)
        
        # Create an empty dictionary to store the data
        data_dict = {}
        
        # Loop through each dataset
        for key in hdf5_file:
            # Select the current dataset
            dataset = hdf5_file[key]
            # Get the data from the dataset
            dataset_data = np.array(dataset)
            # Add the data to the dictionary with the dataset name as the key
            data_dict[key] = dataset_data

        # Close the HDF5 file
        hdf5_file._close()

        return data_dict

    def calculate_P90(self, partial, pressure, Reversed_Order=False):
        # In a couple of cases (Altzomoni and the Azores, specifically) the order of values is mysteriously reversed.
        if Reversed_Order:
            partial = partial[::-1]
            pressure = pressure[::-1]
            
        target = 0.9*sum(partial)
        # Do this the dumb way: start right before the end and keep going until you get bigger than the target
        end = partial.size
        ii = end - 2
        while(sum(partial[ii:end]) < target):
            ii-=1
            # Once you have a value bigger than the target, linearly interpolate between the two pressure values
        # m = (y2-y1)/(x2-x1)
        m = (pressure[ii+1] - pressure[ii])/(sum(partial[(ii+1):end]) - sum(partial[ii:end]))
        # y3 = y1 + m(x3-x1)
        interpolated_P90 = pressure[ii] + m*(target-sum(partial[ii:end]))
        return interpolated_P90

    def P90_search(self, partial, pressure, Reversed_Order=False):

        # In a couple of cases (Altzomoni and the Azores, specifically) the order of values is mysteriously reversed.
        if Reversed_Order:
            partial = partial[::-1]
            pressure = pressure[::-1]

        target = 0.9*sum(partial)
        # Do this the dumb way: start right before the end and keep going until you get bigger than the target
        end = partial.size
        ii = end - 2
        while(sum(partial[ii:end]) < target):
            ii-=1

        # Once you have a value bigger than the target, linearly interpolate between the two pressure values
        # m = (y2-y1)/(x2-x1)
        m = (pressure[ii+1] - pressure[ii])/(sum(partial[(ii+1):end]) - sum(partial[ii:end]))
        # y3 = y1 + m(x3-x1)
        interpolated_P90 = pressure[ii] + m*(target-sum(partial[ii:end]))
        return interpolated_P90

    def add_P90(self, dictionary, Reversed_P):
        #print("Reversed_P value = {}".format(Reversed_P))
        P90_r = np.full(len(dictionary['date']), np.nan)
        P90_ap = np.full(len(dictionary['date']), np.nan)
        for ii in range(len(dictionary['date'])):
            pressure = dictionary['PRESSURE_INDEPENDENT'][ii,]
            column = dictionary['C2H6.COLUMN_ABSORPTION.SOLAR'][ii]
            partial_r = dictionary['C2H6.COLUMN.PARTIAL_ABSORPTION.SOLAR'][ii,]
            partial_ap = dictionary['C2H6.COLUMN.PARTIAL_ABSORPTION.SOLAR_APRIORI'][ii,]
    
            P90_r[ii] = self.P90_search(partial_r, pressure, Reversed_P)
            P90_ap[ii] = self.P90_search(partial_ap, pressure, Reversed_P)
    
        dictionary['P90'] = P90_r
        dictionary['P90_APRIORI'] = P90_ap
        
        return dictionary

    def return_multiple_neighbors(self, latitudes, longitudes, locations, dlat_grid, dlon_grid):
        # Extract the 9 grid cells +/- 1 grid cell away from each original target.
        # This routine also needs the grid-spacing in lat and lon space (dlat_grid, dlon_grid)
        # and a list to append the values to (to_save)
        neighbor_targets=[]
        lat_boundaries = [-90 + dlat_grid, 90 - dlat_grid]
        lon_boundaries = [-180 + dlon_grid, 179.375 - dlon_grid]
        # Check if latitudes are outside the polar grid cells
        if any(lat <= lat_boundaries[0] or lat >= lat_boundaries[1] for lat in latitudes):
            print("Woah, can't return neighbors of the polar grid cells")
            return None


        for lat, lon, location in zip(latitudes, longitudes, locations):
            if lon <= lon_boundaries[0]:
                targets = [(lat + dlat_grid, lon, location),
                           (lat + dlat_grid, lon + dlon_grid, location),
                           (lat + dlat_grid, lon_boundaries[1] - (dlon_grid / 2), location),
                           (lat, lon, location),
                           (lat, lon + dlon_grid, location),
                           (lat, lon_boundaries[1] - (dlon_grid / 2), location),
                           (lat - dlat_grid, lon, location),
                           (lat - dlat_grid, lon + dlon_grid, location),
                           (lat - dlat_grid, lon_boundaries[1] - (dlon_grid / 2), location)]
            elif lon >= lon_boundaries[1]:
                targets = [(lat + dlat_grid, lon, location),
                           (lat + dlat_grid, lon_boundaries[0] + (dlon_grid / 2), location),
                           (lat + dlat_grid, lon - dlon_grid, location),
                           (lat, lon, location),
                           (lat, lon_boundaries[0] - (dlon_grid / 2), location),
                           (lat, lon - dlon_grid, location),
                           (lat - dlat_grid, lon, location),
                           (lat - dlat_grid, lon_boundaries[0] - (dlon_grid / 2), location),
                           (lat - dlat_grid, lon - dlon_grid, location)]
            else:
                targets = [(lat + dlat, lon + dlon, location)
                           for dlat, dlon in [(0, dlon_grid), (0, 0), (0, -dlon_grid), (dlat_grid, dlon_grid), (dlat_grid, 0),
                                              (dlat_grid, -dlon_grid), (-dlat_grid, dlon_grid), (-dlat_grid, 0), (-dlat_grid, -dlon_grid)]]

                neighbor_targets.extend(targets)

        # Now iterate through this and append the values.
        # I don't love that we're not returning anything and instead modifying the input but
        # that may have to be how this works for now
        neighbor_df = pd.DataFrame(neighbor_targets, columns=['Latitudes', 'Longitudes', 'Locations'])

        return neighbor_df

    def select_neighbors(self, df, latitudes, longitudes, dlat_grid, dlon_grid, to_save):
        # for a given xarray dataset df and a set of paired lat-lon indicies (latitudes, longitudes)
        neighbor_df = self.return_multiple_neighbors(latitudes, longitudes, dlat_grid, dlon_grid)
        Location_Selected = df.sel(lat = neighbor_df['Latitudes'].values, lon = neighbor_df['Longitudes'].values,
                                   method = 'nearest')
        to_save.append(Location_Selected)

    def select_one_locations_neighbors(self, df, lat, lon, dlat_grid, dlon_grid):
        # Same thing as for select_neighbors but it works for a single location instead of multiple.
        # I hate making two of these but am not good enough at python to make the general version.
        neighbor_targets=[]
        lat_boundaries = [-90 + dlat_grid, 90 - dlat_grid]
        lon_boundaries = [-180 + dlon_grid, 179.375 - dlon_grid]

        # The routine breaks flat out if you use it at the poles

        if(lat <= lat_boundaries[0] or lat >= lat_boundaries[1]):
            print("Woah, can't return neighbors of the polar gridcells")
            return null

        # Handle the edge cases of the overlaps
        if(lon <= lon_boundaries[0]):
            neighbor_targets.append([(lat + dlat_grid, lon)])
            neighbor_targets.append([(lat + dlat_grid, lon + dlon_grid)])
            neighbor_targets.append([(lat + dlat_grid, lon_boundaries[1]-(dlon_grid/2))])
            neighbor_targets.append([(lat, lon)])
            neighbor_targets.append([(lat, lon + dlon_grid)])
            neighbor_targets.append([(lat, lon_boundaries[1]-(dlon_grid/2))])
            neighbor_targets.append([(lat - dlat_grid, lon)])
            neighbor_targets.append([(lat - dlat_grid, lon + dlon_grid)])
            neighbor_targets.append([(lat - dlat_grid, lon_boundaries[1]-(dlon_grid/2))])

        elif(lon >= lon_boundaries[1]):
            neighbor_targets.append([(lat + dlat_grid, lon)])
            neighbor_targets.append([(lat + dlat_grid, lon_boundaries[0]+(dlon_grid/2))])
            neighbor_targets.append([(lat + dlat_grid, lon - dlon_grid)])
            neighbor_targets.append([(lat, lon)])
            neighbor_targets.append([(lat, lon_boundaries[0]-(dlon_grid/2))])
            neighbor_targets.append([(lat, lon - dlon_grid)])
            neighbor_targets.append([(lat - dlat_grid, lon)])
            neighbor_targets.append([(lat - dlat_grid, lon_boundaries[0]-(dlon_grid/2))])
            neighbor_targets.append([(lat - dlat_grid, lon - dlon_grid)])

        # normal behavior
        else:
            neighbor_targets += [
                (lat + dlat, lon + dlon)
                for dlat, dlon in [(0, dlon_grid), (0, 0), (0, -dlon_grid), (dlat_grid, dlon_grid), (dlat_grid, 0),
                                   (dlat_grid, -dlon_grid), (-dlat_grid, dlon_grid), (-dlat_grid, 0), (-dlat_grid, -dlon_grid)]
            ]
        # Now iterate through this and append the values.
        # I don't love that we're not returning anything and instead modifying the input but 
        # that may have to be how this works for now
        neighbor_df = pd.DataFrame(neighbor_targets, columns=['Latitudes', 'Longitudes'])
        #print(neighbor_df)
        Location_Selected = df.sel(lat = neighbor_df['Latitudes'].values, lon = neighbor_df['Longitudes'].values,
                                   method = 'nearest')
        #to_save.append(Location_Selected)
        return Location_Selected


    def check_and_remove_files(file_names):
        # This subroutine will check the list of files expected and remove non-existent filenames
        # This is necessary for the parallel open function, since otherwise it all fails.
        existing_files = []

        for file_name in file_names:
            if os.path.isfile(file_name):
                existing_files.append(file_name)

        return existing_files

    def Process_NDACC(self, data_dict, datetime_name, variables_to_save, Reversed_P, P90_Flag=False):
        """
        This method takes a data_dict generated with the Read_HDF4 method
        It first converts the MDJ2k to a date and time formate
        Then it formats the output as a pandas dataframe with columns for date, time, and specified variables
        """
        
        # First, convert the date
        data_dict['date'] = []
        for MDJ22K in data_dict['DATETIME']: data_dict['date'].append(MDJ22K_Converter.mj2k2gregorianday(MDJ22K))

        if P90_Flag:
            data_dict = self.add_P90(data_dict, Reversed_P) 

        variables_to_save.insert(0, "date")
        # Now select the relevant key pairs:
        selected_values = {key:data_dict[key] for key in variables_to_save}
        
        # Gotta fix the type of 'date' manually, not certain why
        selected_values['date'] = np.array(selected_values['date'])

        # Now convert to the PANDAS dataframe
        selected_df = pd.DataFrame(selected_values)
        selected_df['date'] = pd.to_datetime(selected_df['date'])
        selected_df['day'] = selected_df['date'].dt.date
        selected_df['time'] = selected_df['date'].dt.time

        return selected_df
        
    def select_value(self, df, overpass_time='13:30'):
        # Sort the DataFrame by 'day' and 'time'
        df = df.sort_values(by=['day', 'time'])

        # Now create the Overpass Time variable for each day
        df['overpass_time'] = pd.to_datetime(df['day'].astype(str)[0] + " " + '13:30')

        # Group the DataFrame by 'day'
        grouped = df.groupby('day')

        # Initialize an empty list to store the selected values
        selected_values = []

        # Iterate over each group
        for day, group in grouped:
            # Check if the group has multiple 'time' entries
            if len(group) > 1:
                # Calculate the time difference from each 'time' to 13:30
                group['Time_Diff'] = abs(group['time'] - pd.to_datetime(overpass_time).time())

                # Sort the group by 'time_Diff' in ascending order
                group = group.sort_values(by='Time_Diff')

                # Select the value nearest to 13:30
                selected_value = group.iloc[0]['Value']
            else:
                # If the group has only one entry, select the value directly
                selected_value = group.iloc[0]['Value']

            # Append the selected value to the list
            selected_values.append(selected_value)

        # Create a new DataFrame with 'day' and the selected values
        result_df = pd.DataFrame({'day': grouped['day'].first(), 'Selected_Value': selected_values})

        return result_df
