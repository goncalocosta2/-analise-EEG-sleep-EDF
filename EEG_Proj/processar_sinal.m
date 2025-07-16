function [X, y] = processar_sinal(ficheiro_psg, ficheiro_hyp, filtro)
    fs = 100;  % Frequência de amostragem (Hz)
    janela_seg = 30;
    tam_janela = fs * janela_seg;

    % 1. Carregar dados EDF
    psg = edfread(ficheiro_psg);
    hyp = edfread(ficheiro_hyp);
    info = edfinfo(ficheiro_hyp);

    % 2. Verificar canal
    canal_nome = 'EEGFpz_Cz';
    eeg = psg.(canal_nome);

    % Corrigir tipo de dado (caso seja cell)
    if iscell(eeg)
        try
            eeg = cell2mat(eeg);
        catch
            X = []; y = [];
            return;
        end
    end
    eeg = double(eeg);  % Garantir double

    % 3. Aplicar filtro
    eeg = conv(eeg, filtro, 'same');

    % 4. Anotações 
    anot = info.Annotations;
    estagios_txt = anot.Annotations;
    onsets = seconds(anot.Onset);
    duracoes = seconds(anot.Duration);

    dic_estagios = containers.Map( ...
        {'Sleep stage W', 'Sleep stage 1', 'Sleep stage 2', ...
         'Sleep stage 3', 'Sleep stage 4', 'Sleep stage R'}, ...
        [1, 2, 3, 4, 4, 6]);  % N3 = estágios 3 e 4

    n_janelas = floor(length(eeg) / tam_janela);
    X = [];
    y = [];

    % 5. Definir bandas e vetores
    n_fft = 1024;
    f = (0:n_fft-1)*(fs/n_fft);
    freqs = f(1:n_fft/2);

    bandas = {
        'Delta', [0.5 4];
        'Teta',  [4 8];
        'Alfa',  [8 13];
        'Beta',  [13 30];
    };

    % 6. Loop por janelas
    for j = 1:n_janelas
        ini = (j-1)*tam_janela + 1;
        fim = j*tam_janela;
        janela = eeg(ini:fim);

        % Timestamp da janela
        tempo_janela = (ini-1)/fs;

        % Verificar anotação válida
        idx_estagio = find(tempo_janela >= onsets & tempo_janela < (onsets + duracoes), 1);
        if isempty(idx_estagio), continue; end

        nome_estagio = estagios_txt{idx_estagio};
        if ~isKey(dic_estagios, nome_estagio), continue; end
        rotulo = dic_estagios(nome_estagio);

        % Normalização local
        janela = (janela - mean(janela)) / (std(janela) + eps);

        % FFT e PSD
        fft_janela = fft(janela, n_fft);
        psd = abs(fft_janela(1:n_fft/2)).^2 / n_fft;

        potencias = zeros(1, length(bandas));
        for b = 1:length(bandas)
            faixa = bandas{b,2};
            idx_freqs = freqs >= faixa(1) & freqs <= faixa(2);
            potencias(b) = sum(psd(idx_freqs));
        end

        total = sum(psd);
        delta_teta = potencias(1) / (potencias(2) + eps);

        % Features adicionais
        media = mean(janela);
        desvio = std(janela);
        energia = sum(janela.^2);
        zcr = sum(abs(diff(janela > 0))) / length(janela);

        % Vetor final de características
        carac = [potencias, total, delta_teta, media, desvio, energia, zcr];

        X = [X; carac];
        y = [y; rotulo];
    end
end
