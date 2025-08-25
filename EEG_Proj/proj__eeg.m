clc;
clear all;
close all;

%% 1. Carregar Dados

PSG='SC4021E0-PSG.edf';
HYP='SC4021EH-Hypnogram.edf';

% Ler o sinal EEG
psg_dados = edfread(PSG); 
eeg_sinal = psg_dados.EEGFpz_Cz;  % Pega o sinal EEG da coluna correta
fs = 100;  % Frequência de amostragem (Hz)

% Ler o hipnograma
hyp_dados = edfread(HYP);
info = edfinfo(HYP);
anotacoes = info.Annotations;
tempo = seconds(anotacoes.Onset);
estagios = anotacoes{:,"Annotations"};

%% 2. Pré-processamento dos Estágios (mapear para janelas de 30s)

% Dicionário dos estágios
dic_estagios = containers.Map(... 
    {'Sleep stage W', 'Sleep stage 1', 'Sleep stage 2', 'Sleep stage 3', ...
    'Sleep stage 4', 'Sleep stage R'}, ...
    {1, 2, 3, 4, 5, 6});

% Duração de cada janela
janela_seg = 30;
duracao_total = seconds(anotacoes.Onset(end)) + seconds(anotacoes.Duration(end));
num_janelas = floor(duracao_total / janela_seg);

% Vetor de estágios 
estagios_cod = zeros(num_janelas, 1);

% Preencher o vetor de estágios janela a janela
for i = 1:height(anotacoes)
    nome_estagio = anotacoes.Annotations{i};
    
    if isKey(dic_estagios, nome_estagio)
        estagio_cod = dic_estagios(nome_estagio);
    else
        continue;
    end

    inicio = seconds(anotacoes.Onset(i));
    duracao = seconds(anotacoes.Duration(i));
    fim = inicio + duracao;
    
    % Determinar janelas que caem nesse intervalo
    idx_inicio = floor(inicio / janela_seg) + 1;
    idx_fim = floor(fim / janela_seg);
    
    % Proteger contra valores fora dos limites
    idx_inicio = max(1, idx_inicio);
    idx_fim = min(num_janelas, idx_fim);

    estagios_cod(idx_inicio:idx_fim) = estagio_cod;
end

%% 3. Filtro Passa-Banda (0.5 Hz - 45 Hz) 
% Passar o formato para double
if iscell(eeg_sinal)
    eeg_sinal = cell2mat(eeg_sinal);
end
eeg_sinal = double(eeg_sinal);  

% Filtro passa-banda entre 0.5 Hz e 45 Hz com fir1

ordem = 500;  % Ordem do filtro (quanto maior, mais seletivo)
frequencia_corte = [0.5 45] / (fs/2);  % Normalizar pela frequência de Nyquist

% Criar filtro FIR passa-banda com janela Hamming
filtro = fir1(ordem, frequencia_corte, 'bandpass', hamming(ordem + 1));

%% 4. Filtragem do sinal EEG

% Aplicar o filtro com convolução 
eeg_sinal_filtrado = filtfilt(filtro, 1, eeg_sinal); 

% Plotar o sinal original e o sinal filtrado
figure;

% Sinal original
subplot(2, 1, 1);
plot(eeg_sinal);
xlabel('Índice de Amostra');
ylabel('Amplitude');
title('Sinal EEG Original');
grid on;

% Sinal filtrado
subplot(2, 1, 2);
plot(eeg_sinal_filtrado);
xlabel('Índice de Amostra');
ylabel('Amplitude');
title('Sinal EEG Filtrado');
grid on;

%% 5. Análise Espectral (Distribuição da Potência)

% Parâmetros da FFT
n_fft = 1024;  % Número de pontos da FFT sem exigir muito esforço do pc
tam_janela = 30 * fs;  % 30 segundos de janela porque cada estagio de sono tem 30 segundos 
n_amostras = length(eeg_sinal_filtrado);
f = (0:n_fft-1)*(fs/n_fft);  % Frequências (Hz)

% Dividir o sinal em janelas e calcular a FFT 
num_janelas = floor(n_amostras / tam_janela);  % Número de janelas possíveis
P_eeg = zeros(num_janelas, n_fft / 2);  % Matriz para armazenar as potências de frequência 

for j = 1:num_janelas
    % Definir o início e o fim da janela
    inicio = (j-1) * tam_janela + 1;
    fim = min(j * tam_janela, n_amostras);
    
    janela = eeg_sinal_filtrado(inicio:fim);
    
    % Transformar o formato da janela em double caso este não esteja
    janela = double(janela);
    
    janela = janela - mean(janela);       % Remover componente média
    janela = janela / std(janela);        % Normalizar a amplitude

    
    % Calcular a FFT e a potência espectral
    fft_sinal = fft(janela, n_fft);
    N = length(janela);
    P_sinal = (1/(fs*N)) * abs(fft_sinal(1:n_fft/2)).^2;
    P_sinal(2:end-1) = 2*P_sinal(2:end-1); % corrigir parte positiva
    
    % Armazenar a potência da janela na matriz P_eeg
    P_eeg(j, :) = P_sinal;
end

%% 6. Visualização: Distribuição da Potência

% Plotar PSD (Distribuição da Potência) para a primeira janela
figure;
plot(f(1:n_fft/2), P_eeg(1, :)); 
xlabel('Frequência (Hz)');
ylabel('Densidade espectral de potência');
title('Distribuição da Potência (PSD) - Primeira Janela (30s)');
grid on;

figure;
plot(f(1:n_fft/2), abs(fft_sinal(1:n_fft/2)));
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('FFT - Magnitude da Primeira Janela (30s)');
grid on;


%% 7. Cálculo da Potência por Banda de Frequência

% Definir as bandas de frequência (em Hz)
bandas_freq = {
    'Delta', [0.5 4];
    'Teta',  [4 8];
    'Alfa',  [8 13];
    'Beta',  [13 30];
    'Gama', [30 40];
};

% Frequências associadas à FFT (metade positiva)
freqs = f(1:n_fft/2);

% Ajustar o vetor de estágios para ter o mesmo número de janelas que P_eeg
estagios_cod = estagios_cod(1:size(P_eeg,1));

% Obter a lista de estágios únicos
estagios_unicos = unique(estagios_cod);


% Inicializar matriz para armazenar a potência média por banda e por estágio
pot_por_banda = zeros(length(estagios_unicos), size(bandas_freq, 1));

% Loop pelos diferentes estágios de sono
for idx_e = 1:length(estagios_unicos)
    estagio = estagios_unicos(idx_e);
    
    % Encontrar as janelas que pertencem a esse estágio
    idx_janelas = find(estagios_cod == estagio);
    
    % Calcular média espectral (PSD) das janelas desse estágio
    media_psd = mean(P_eeg(idx_janelas, :), 1);
    
    % Calcular a potência em cada banda de frequência
    for idx_b = 1:size(bandas_freq, 1)
        faixa = bandas_freq{idx_b, 2};
        indices_banda = freqs >= faixa(1) & freqs <= faixa(2);
        df = fs/n_fft;
        pot_banda = sum(media_psd(indices_banda)) * df;
        pot_por_banda(idx_e, idx_b) = pot_banda;
    end
end

% Normalizar a potência por banda (em % da potência total de cada banda)
pot_por_banda = pot_por_banda ./ sum(pot_por_banda, 1);
pot_por_banda = pot_por_banda * 100;

%% 8. Visualização da Potência Relativa por Banda e Estágio

% Gráfico de barras
figure;
bar(pot_por_banda, 'grouped');
xlabel('Estágio do Sono');
ylabel('Potência Relativa (%)');
title('Distribuição da Potência por Banda em Cada Estágio de Sono');
xticklabels({'W', 'N1', 'N2', 'N3', 'N4', 'REM'});  
legend({'Delta', 'Teta', 'Alfa', 'Beta', 'Gama'});
grid on;


%% 9. Potencia Média por Banda e Estágio 

% Nomes das bandas
nomes_bandas= {'Delta', 'Teta', 'Alfa', 'Beta', 'Gama'};
num_bandas = length(nomes_bandas);

% Vetores para características e rótulos
caracteristicas = [];
rotulos = [];

for j = 1:num_janelas
    potencias_banda = zeros(1, num_bandas);

    for b = 1:num_bandas
        faixa = bandas_freq{b, 2};
        indices_freq = freqs >= faixa(1) & freqs <= faixa(2);
        potencias_banda(b) = sum(P_eeg(j, indices_freq));
    end

    caracteristicas = [caracteristicas; potencias_banda];

    % Associar o estágio à janela j
    if j <= length(estagios_cod)
        rotulos = [rotulos; estagios_cod(j)];
    else
        rotulos = [rotulos; 0];  % Caso não tenha rótulo
    end
end

% Estágios ordenados
num_estagios = length(estagios_unicos);

% Inicializar matriz de médias 
media_por_estagio = zeros(num_estagios, num_bandas);

% Calcular a média das potências por banda para cada estágio
for i = 1:num_estagios
    idx = rotulos == estagios_unicos(i);
    media_por_estagio(i, :) = mean(caracteristicas(idx, 1:num_bandas), 1);
end

%% 10 Grafico da Potencia média 

% Gráfico de barras com escala logarítmica
figure;
bar(nomes_bandas, media_por_estagio');
xlabel('Bandas de Frequência');
ylabel('Potência Média');
title('Potência Média por Banda para Cada Estágio do Sono');
set(gca, 'YScale', 'log');  % Escala logarítmica para melhor visualização
legend({'Sleep stage W', 'Sleep stage 1', 'Sleep stage 2', ...
        'Sleep stage 3', 'Sleep stage 4', 'Sleep stage R'}, 'Location', 'northwest');
grid on;

%% 11. Extração de Características para Classificação do jogo

% Estágio aleatório entre os disponíveis
estagio_escolhido = estagios_unicos(randi(length(estagios_unicos)));

% Janelas com esse estágio
idx_estagio = find(rotulos == estagio_escolhido);

% Janela aleatória entre as que têm esse estágio
idx_central = idx_estagio(randi(length(idx_estagio)));

% Procurar para trás
inicio = idx_central;
while inicio > 1 && rotulos(inicio - 1) == estagio_escolhido
    inicio = inicio - 1;
end

% Procurar para frente
fim = idx_central;
while fim < length(rotulos) && rotulos(fim + 1) == estagio_escolhido
    fim = fim + 1;
end

% Janelas todas
bloco_indices = inicio:fim;

fprintf('Bloco de janelas: %d até %d (%d janelas)\n', ...
     inicio, fim, length(bloco_indices));

% Potência média por banda nesse bloco
potencias_banda_bloco = zeros(1, size(bandas_freq,1));
for b = 1:size(bandas_freq,1)
    faixa = bandas_freq{b, 2};
    indices_freq = freqs >= faixa(1) & freqs <= faixa(2);
    
    % Potência média da banda no bloco
    potencias_banda_bloco(b) = mean(sum(P_eeg(bloco_indices, indices_freq), 2));
end


disp('Potência média por banda no bloco:');
for b = 1:size(bandas_freq,1)
    fprintf('%s: %.2f\n', bandas_freq{b,1}, potencias_banda_bloco(b));
end

figure;
bar(categorical(bandas_freq(:,1)), potencias_banda_bloco);
title(sprintf('Potência Média por Banda '));
ylabel('Potência');
grid on;


%12 Jogo: Adivinha o Estágio de Sono

% Mostrar opções ao utilizador
fprintf('\nCom base nestas potências médias por banda...\n');
fprintf('Qual achas que é o estágio de sono?\n');
fprintf('1 - W (vigília)\n2 - N1\n3 - N2\n4 - N3\n5 - N4\n6 - REM\n');

% Loop até acertar
acertou = false;
tentativas = 0;

while ~acertou
    resposta = input('Insere o número correspondente ao estágio: ');
    tentativas = tentativas + 1;

    if resposta == estagio_escolhido
        acertou = true;
        fprintf('\n Muito bem! Acertaste o estágio na %dª tentativa!\n', tentativas);
    else
        fprintf(' Errado. Tenta outra vez!\n');
    end
end

% Descrição de cada estágio
descricoes = {
    'W',  'Acertaste o estágio 1, que corresponde à vigília (W). Neste estado, o EEG mostra atividade de alta frequência e baixa amplitude, principalmente nas bandas alfa e beta. O cérebro está alerta e os olhos podem estar abertos ou a piscar.';
    'N1', 'Acertaste o estágio 2, correspondente ao N1. Este é o primeiro estágio do sono leve. Caracteriza-se por uma redução da atividade alfa e aumento na banda teta. O EEG mostra transições suaves e o corpo começa a relaxar.';
    'N2', 'Acertaste o estágio 3, correspondente ao N2. Este estágio é marcado por atividade predominante na banda teta e pela presença de fusos do sono, sinais típicos neste estágio intermediário do sono leve.';
    'N3', 'Acertaste o estágio 4, correspondente ao N3. É um estágio de sono profundo, com predominância de ondas lentas (delta). O EEG mostra alta amplitude e baixa frequência. Essencial para a recuperação física e consolidação da memória.';
    'N4', 'Acertaste o estágio 5, correspondente ao N4. É um estágio mais profundo que N3, com uma predominância ainda maior da banda delta. Muitas vezes é agrupado com o N3 nas classificações atuais.';
    'REM','Acertaste o estágio 6, correspondente ao sono REM. Apesar de ser um estágio de sono profundo, o EEG apresenta padrões semelhantes à vigília, com atividade mista nas bandas beta e teta. É o estágio em que ocorrem os sonhos mais intensos.';
};

% Mostrar descrição correspondente
%nome_estagio = descricoes{estagio_escolhido, 1};
texto_descricao = descricoes{estagio_escolhido, 2};

fprintf('%s\n', texto_descricao);

%% 13.Experiencia de um Classificador com TreeBagger

% Lista de ficheiros PSG e Hypnogram para treino
ficheiros = {
    'SC4001E0-PSG.edf', 'SC4001EC-Hypnogram.edf';
    'SC4011E0-PSG.edf', 'SC4011EH-Hypnogram.edf';
    'SC4021E0-PSG.edf', 'SC4021EH-Hypnogram.edf';
    'SC4031E0-PSG.edf', 'SC4031EC-Hypnogram.edf';
    'SC4051E0-PSG.edf', 'SC4051EC-Hypnogram.edf';
    'SC4061E0-PSG.edf', 'SC4061EC-Hypnogram.edf';
};

% Dados de treino
X_treino = [];
y_treino = [];

for k = 1:size(ficheiros,1)
    [X_tmp, y_tmp] = processar_sinal(ficheiros{k,1}, ficheiros{k,2}, filtro);
    X_treino = [X_treino; X_tmp];
    y_treino = [y_treino; y_tmp];
end

% Normalização
X_treino = normalize(X_treino);

% PCA (redução para 5 componentes)
[coeff, X_treino_pca, ~, ~, explained] = pca(X_treino);
X_treino_pca = X_treino_pca(:,1:5);

% Treinar TreeBagger
modelo = TreeBagger(100, X_treino_pca, y_treino, 'OOBPrediction','On', 'Method','classification');

%% 14. Avaliação em Múltiplos Ficheiros de Teste

ficheiros_teste = {
    'SC4022E0-PSG.edf', 'SC4022EJ-Hypnogram.edf';
    'SC4002E0-PSG.edf', 'SC4002EC-Hypnogram.edf';
};

for i = 1:size(ficheiros_teste,1)
    PSG = ficheiros_teste{i,1};
    HYP = ficheiros_teste{i,2};

    % Processar sinal
    [X_teste, y_teste] = processar_sinal(PSG, HYP, filtro);
    X_teste = normalize(X_teste);
    X_teste_pca = X_teste * coeff(:,1:5);  % Aplicar o mesmo PCA

    % Previsões
    y_pred = str2double(predict(modelo, X_teste_pca));

    % Precisão
    acc = sum(y_pred == y_teste) / length(y_teste) * 100;
    fprintf('Precisão para %s: %.2f%%\n', PSG, acc);

    % Matriz de confusão
    figure;
    confusionchart(y_teste, y_pred, ...
        'RowSummary','row-normalized', ...
        'ColumnSummary','column-normalized');
    title(['Matriz de Confusão - ', PSG]);
end
