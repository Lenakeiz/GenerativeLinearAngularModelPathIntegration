
HOMCI_model = load('Data/Delong_Model_HOMCI.mat');
HOMCI_data_linear = load('Data/Delong_Data_HOMCI_Linear.mat');
HOMCI_data_angular = load('Data/Delong_Data_HOMCI_Angular.mat');



HOMCI.ratings = [HOMCI_model.sample.ratings;
            HOMCI_data_linear.sample.ratings;
            HOMCI_data_angular.sample.ratings];
HOMCI.spsizes = HOMCI_model.sample.spsizes;

%save
save('Data/HOMCI.mat', '-struct', 'HOMCI')



MCI_model = load('Data/Delong_Model_MCI.mat');
MCI_data_linear = load('Data/Delong_Data_MCI_Linear.mat');
MCI_data_angular = load('Data/Delong_Data_MCI_Angular.mat');



MCI.ratings = [MCI_model.sample.ratings;
            MCI_data_linear.sample.ratings;
            MCI_data_angular.sample.ratings];
MCI.spsizes = MCI_model.sample.spsizes;

%save
save('Data/MCI.mat', '-struct', 'MCI')