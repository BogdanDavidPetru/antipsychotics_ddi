import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

drug_interaction_df = pd.read_csv('drug_interaction_train_test_data_not_valid_antipsychotics.csv')

print("----------Initial Data Frame info-------------")
drug_interaction_df = drug_interaction_df.iloc[:, 1:]  # remove first column which represents index
print(drug_interaction_df.info())

print("----------Input Data Frame info-------------")
input_drug_interaction_df = drug_interaction_df.iloc[:, :-1]
print(input_drug_interaction_df.info())

print("----------Output Data Frame info-------------")
output_drug_interaction_df = drug_interaction_df.iloc[:, -1:]
print(output_drug_interaction_df.info())

X_train, X_test, y_train, y_test = train_test_split(input_drug_interaction_df, output_drug_interaction_df,
                                                    test_size=0.3, random_state=42)

X_tr_arr = X_train
X_ts_arr = X_test
y_tr_arr = y_train.values
y_ts_arr = y_test.values

print('Input Shape', X_tr_arr.shape)
print('Output Shape', X_test.shape)
print('Output train shape: ', y_tr_arr.shape)
print('Output test shape: ', y_ts_arr.shape)

#train

regression_model = LogisticRegression(solver='lbfgs', max_iter=1000)

regression_model.fit(X_tr_arr, y_tr_arr.ravel())

pred = regression_model.predict(X_ts_arr)

print('Accuracy from sk-learn: {0}'.format(regression_model.score(X_ts_arr, y_ts_arr.ravel())))

