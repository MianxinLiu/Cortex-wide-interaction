import pandas as pd
import numpy as np
import scipy.io as scio

import torch
import torch.nn as nn
from torch.utils.data import Dataset, TensorDataset
import MyModels
from sklearn import metrics
from sklearn.model_selection import KFold
from random import shuffle
import mat73
import glob

for start in [-25]:
    permu_results = []

    X = mat73.loadmat('./Data/pre_vol.mat')
    X = X['Vol']

    temp = scio.loadmat('./Python/shuffled_index_all')
    index = temp['index'][0]

    end = 0
    #end = start+25
    seglen = end - start
    S = 500 + start
    E = 500 + end
    X = X[:, S:E, index]
    X = np.transpose(X, [2, 1, 0])
    X = torch.from_numpy(X)
    X = torch.tensor(X, dtype=torch.float32)

    for r in range(100):
        Y = scio.loadmat('./Data/label_hie_all.mat')
        Y = Y['label']
        shuffle(Y[0, :])
        # Y=Y[:,index]
        Y = np.squeeze(Y)
        Y = torch.from_numpy(Y)
        Y = Y - 1
        Y = torch.tensor(Y, dtype=torch.float32)

        qual_all = []
        for cv in [1, 2, 3, 4, 5]:
            temp = scio.loadmat('./Python//kfold/all/shuffled_index_cv' + str(cv))
            train_idx = temp['train_idx'][0]
            test_idx = temp['test_idx'][0]
            print("cv:", cv)
            dataset_train = TensorDataset(X[train_idx, :, :], Y[train_idx])
            dataset_test = TensorDataset(X[test_idx, :, :], Y[test_idx])

            train_loader = torch.utils.data.DataLoader(dataset_train, batch_size=30, shuffle=True, drop_last=True)
            test_loader = torch.utils.data.DataLoader(dataset_test, batch_size=30)

            ttmp = Y[train_idx]
            ratio1 = ttmp[ttmp == 1].shape[0] / ttmp[ttmp == 0].shape[0]
            ratio2 = ttmp[ttmp == 2].shape[0] / ttmp[ttmp == 0].shape[0]
            weight = torch.cuda.FloatTensor([1, 1 / ratio1, 1 / ratio2])
            # weight = torch.cuda.FloatTensor([1, 1.5])
            loss_func = nn.CrossEntropyLoss(weight)  # the target label is not one-hotted

            lr = 0.001
            EPOCH = 50

            test_auc = []
            train_los = []
            test_los = []
            train_auc = []
            sen = []
            spe = []
            # epoch_of_lr_decrease = 20

            model = MyModels.RNN(input_size=1236, hidden_size=100, num_layers=1, tsize=seglen, classnum=3)
            model.cuda()
            optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-2)
            # epoch_of_lr_decrease = 20
            auc_baseline = 0.0
            for epoch in range(EPOCH):
                for step, (b_x, b_y) in enumerate(train_loader):  # gives batch data
                    model.train()
                    b_x = b_x.cuda()
                    b_y = b_y.cuda()
                    b_y = b_y.long()
                    output = model(b_x)  # rnn output

                    loss = loss_func(output, b_y)  # cross entropy loss
                    optimizer.zero_grad()  # clear gradients for this training step
                    loss.backward()  # backpropagation, compute gradients
                    optimizer.step()  # apply gradients

                    predicted = torch.max(output.data, 1)[1]
                    correct = (predicted == b_y).sum()
                    tr_accuracy = float(correct) / float(b_x.shape[0])
                    train_auc.append(tr_accuracy)
                    # print('[Epoch %d, Batch %5d] loss: %.3f' %
                    #      (epoch + 1, step + 1, loss))
                    # print('|train diag loss:', loss.data.item(), '|train accuracy:', tr_accuracy
                    #      )
                    if epoch >= 40 and tr_accuracy > 0.8:
                        # lr=0.001
                        # optimizer = torch.optim.Adam(gclstm.parameters(), lr=lr, weight_decay=1e-2)
                        predicted_all = []
                        pre_score_all = []
                        test_y_all = []
                        model.eval()
                        with torch.no_grad():
                            for i, (test_x, test_y) in enumerate(test_loader):
                                test_y = test_y.long()
                                test_y = test_y.cuda()
                                test_x = test_x.cuda()
                                test_output = model(test_x)
                                loss = loss_func(test_output, test_y)
                                # print('[Epoch %d, Batch %5d] valid loss: %.3f' %
                                #      (epoch + 1, step + 1, loss))
                                # test_loss = loss_func(test_output, test_y)
                                predicted = torch.max(test_output.data, 1)[1]
                                correct = (predicted == test_y).sum()
                                accuracy = float(correct) / float(predicted.shape[0])
                                test_y = test_y.cpu()
                                predicted = predicted.cpu()
                                predicted_all = predicted_all + predicted.tolist()
                                pre_score_all = pre_score_all + test_output[:, 1].tolist()
                                test_y_all = test_y_all + test_y.tolist()

                        correct = (np.array(predicted_all) == np.array(test_y_all)).sum()
                        accuracy = float(correct) / float(len(test_y_all))
                        test_auc.append(accuracy)
                        sen_1 = metrics.recall_score(test_y_all, predicted_all, labels=[0], average='macro')
                        sen_2 = metrics.recall_score(test_y_all, predicted_all, labels=[1], average='macro')
                        sen_3 = metrics.recall_score(test_y_all, predicted_all, labels=[2], average='macro')
                        pre_1 = metrics.precision_score(test_y_all, predicted_all, labels=[0], average='macro')
                        pre_2 = metrics.precision_score(test_y_all, predicted_all, labels=[1], average='macro')
                        pre_3 = metrics.precision_score(test_y_all, predicted_all, labels=[2], average='macro')
                        F1 = metrics.f1_score(test_y_all, predicted_all, average='macro')

                        if accuracy >= auc_baseline and sen_1 >= 0.0 and sen_2 >= 0.0 and sen_3 >= 0.0:
                            auc_baseline = accuracy
                            torch.save(model.state_dict(), './Python/model_v/temp25_cv' + str(cv) + '.pth')

        qual_all = []
        predicted_all = []
        pre_score_all = []
        test_y_all = []
        for cv in [1, 2, 3, 4, 5]:
            temp = scio.loadmat('./Python/kfold/all/shuffled_index_cv' + str(cv))
            test_idx = temp['test_idx'][0]
            dataset_test = TensorDataset(X[test_idx, :, :], Y[test_idx])
            test_loader = torch.utils.data.DataLoader(dataset_test, batch_size=30)

            model = MyModels.RNN(input_size=1236, hidden_size=100, num_layers=1, tsize=seglen, classnum=3)
            model.load_state_dict(torch.load(
                './Python/model_v/temp25_cv' + str(cv) + '.pth'))
            model.cuda()
            model.eval()
            with torch.no_grad():
                for i, (test_x, test_y) in enumerate(test_loader):
                    test_y = test_y.long()
                    test_y = test_y.cuda()
                    test_x = test_x.cuda()
                    test_output = model(test_x)
                    predicted = torch.max(test_output.data, 1)[1]
                    correct = (predicted == test_y).sum()
                    accuracy = float(correct) / float(predicted.shape[0])
                    test_y = test_y.cpu()
                    predicted = predicted.cpu()
                    pre_score_all = pre_score_all + test_output[:, 1].tolist()
                    predicted_all = predicted_all + predicted.tolist()
                    test_y_all = test_y_all + test_y.tolist()

        correct = (np.array(predicted_all) == np.array(test_y_all)).sum()
        correct = (np.array(predicted_all) == np.array(test_y_all)).sum()
        accuracy = float(correct) / float(len(test_y_all))
        sen_1 = metrics.recall_score(test_y_all, predicted_all, labels=[0], average='macro')
        sen_2 = metrics.recall_score(test_y_all, predicted_all, labels=[1], average='macro')
        sen_3 = metrics.recall_score(test_y_all, predicted_all, labels=[2], average='macro')
        pre_1 = metrics.precision_score(test_y_all, predicted_all, labels=[0], average='macro')
        pre_2 = metrics.precision_score(test_y_all, predicted_all, labels=[1], average='macro')
        pre_3 = metrics.precision_score(test_y_all, predicted_all, labels=[2], average='macro')
        F1 = metrics.f1_score(test_y_all, predicted_all, average='macro')
        print('realization:', r,
              '|test accuracy:', accuracy,
              '|test acc_1:', sen_1,
              '|test acc_2:', sen_2,
              '|test acc_3:', sen_3,
              '|test pre1:', pre_1,
              '|test pre2:', pre_2,
              '|test pre3:', pre_3,
              '|F1:', F1
              )
        permu_results.append([accuracy, sen_1, sen_2, sen_3, pre_1, pre_2, pre_3, F1])

    print(np.mean(permu_results, axis=0))
    print(np.std(permu_results, axis=0))

    Y = scio.loadmat('./Data/label_hie_all.mat')
    Y = Y['label']
    Y = Y[index, :]
    Y = np.squeeze(Y)
    Y = torch.from_numpy(Y)
    Y = Y - 1
    Y = torch.tensor(Y, dtype=torch.float32)

    qual_all = []
    predicted_all = []
    pre_score_all = []
    test_y_all = []
    for cv in [1, 2, 3, 4, 5]:
        temp = scio.loadmat('./Python/kfold/all/shuffled_index_cv' + str(cv))
        test_idx = temp['test_idx'][0]
        dataset_test = TensorDataset(X[test_idx, :, :], Y[test_idx])
        test_loader = torch.utils.data.DataLoader(dataset_test, batch_size=30)

        model = MyModels.RNN(input_size=1236, hidden_size=100, num_layers=1, tsize=seglen, classnum=3)
        model.load_state_dict(torch.load(
            './Python/model_v/model_all_cv' + str(cv) + '_' + str(start) + '_' + str(end) + '.pth'))
        model.cuda()
        model.eval()
        with torch.no_grad():
            for i, (test_x, test_y) in enumerate(test_loader):
                test_y = test_y.long()
                test_y = test_y.cuda()
                test_x = test_x.cuda()
                test_output = model(test_x)
                predicted = torch.max(test_output.data, 1)[1]
                correct = (predicted == test_y).sum()
                accuracy = float(correct) / float(predicted.shape[0])
                test_y = test_y.cpu()
                predicted = predicted.cpu()
                pre_score_all = pre_score_all + test_output[:, 1].tolist()
                predicted_all = predicted_all + predicted.tolist()
                test_y_all = test_y_all + test_y.tolist()

    correct = (np.array(predicted_all) == np.array(test_y_all)).sum()
    correct = (np.array(predicted_all) == np.array(test_y_all)).sum()
    accuracy = float(correct) / float(len(test_y_all))
    sen_1 = metrics.recall_score(test_y_all, predicted_all, labels=[0], average='macro')
    sen_2 = metrics.recall_score(test_y_all, predicted_all, labels=[1], average='macro')
    sen_3 = metrics.recall_score(test_y_all, predicted_all, labels=[2], average='macro')
    pre_1 = metrics.precision_score(test_y_all, predicted_all, labels=[0], average='macro')
    pre_2 = metrics.precision_score(test_y_all, predicted_all, labels=[1], average='macro')
    pre_3 = metrics.precision_score(test_y_all, predicted_all, labels=[2], average='macro')
    F1 = metrics.f1_score(test_y_all, predicted_all, average='macro')
    qual_all.append([accuracy, sen_1, sen_2, sen_3, pre_1, pre_2, pre_3, F1])

    p_values = np.zeros([8])
    for metric in range(8):
        count = 0
        for r in range(100):
            if permu_results[r][metric] >= qual_all[0][metric]:
                count = count + 1
        p_values[metric] = count / 100
    print(p_values)

    scio.savemat('./Python/permutation/permu'
    +'_'+str(start)+'_'+ str(end) +'_Vol.mat', {'permu_results': permu_results, 'p_values': p_values})